#.(require :alexandria)
(defpackage :kdtree (:use :cl :alexandria))
(in-package :kdtree)
(declaim (optimize (speed 2) (debug 3) (safety 3)))

(defconstant +k+ 2)
(defconstant +big+ 12)
;;(defconstant +max-component+ 2000)
(deftype component ()
  `(double-float))

(deftype vec ()
  `(simple-array component (,+k+)))

(declaim (inline v v- v+ v* v/ v. norm norm2))

(defun v (&optional (x #.(coerce 0 'component))
	  (y #.(coerce 0 'component)) (z #.(coerce 0 'component)))
  (declare (component x y)
	   (values vec &optional))
  (make-array +k+ :element-type 'component
	      :initial-contents (subseq (list x y z) 0 +k+)))

(defmacro with-arrays (arrays &body body)
  "Provides a corresponding accessor for each array as a local macro,
so that (ARRAY ...) corresponds to (AREF ARRAY ...)."
  `(macrolet ,(mapcar (lambda (array)
			`(,array (&rest indices) `(aref ,',array ,@indices)))
		      arrays)
     ,@body))

(defmacro do-vec ((i input-arrays &optional result-array)
		  &body body)
  (let ((args (if result-array
		  `(,@input-arrays ,result-array)
		  `(,@input-arrays))))
   `(let ,(when result-array
		`((,result-array (v))))
      (declare (vec ,@args))
      (with-arrays (,@args)
	(dotimes (,i +k+)
	  ,@body)
	,result-array))))
#+nil
(do-vec (i (a b) c) (setf (c i) (+ (a i) (b i))))
#+nil
(do-vec (i () c) (setf (c i) (+ (a i) (b i))))
#+nil
(do-vec (i ()) (setf (c i) (+ (a i) (b i))))
#+nil
(do-vec (i (b) ) (setf (c i) (+ (a i) (b i))))

(defun v+ (a b)
  (declare (vec a b)
	   (values vec &optional))
  (do-vec (i (a b) c)
    (setf (c i) (+ (a i) (b i)))))

(defun v- (a b)
  (declare (vec a b)
	   (values vec &optional))
  (do-vec (i (a b) c)
    (setf (c i) (- (a i) (b i)))))

(defun v* (a s)
  (declare (vec a)
	   (component s)
	   (values vec &optional))
  (do-vec (i (a) c)
    (setf (c i) (* s (a i)))))

(defun v/ (a s)
  (declare (vec a)
	   (component s)
	   (values vec &optional))
  ;; if we have defined the component type as some kind of float we
  ;; can just multiply with the inverse. otherwise we want a fixnum
  ;; result
  (if #.(sb-kernel::subtypep 'component 'float)
      (v* a (/ s))
      (do-vec (i (a) c)
	(setf (c i) (floor (a i) s)))))
#+nil
(v/ (v 5 2) 3)

(defun v. (a b)
  (declare (vec a b)
	   (values component &optional))
  (let ((sum #.(coerce 0 'component)))
    (do-vec (i (a b))
      (incf sum (* (a i) (b i))))
    sum))
#+nil
(v. (v 1 3) (v 2 3))

(defun norm2 (a)
  (declare (vec a)
	   (values component &optional))
  (v. a a))

(defun norm (a)
  (declare (vec a)
	   (values double-float &optional))
  (sqrt (* 1d0 (norm2 a))))

(defun nearest-naive (points point)
  (declare ((simple-array vec 1) points)
	   (vec point)
	   (values fixnum &optional))
  (with-arrays (points)
    (let* ((nearest 0)
	   (nearest-dist2 (norm2 (v- (points nearest) point))))
      (loop for i from 1 below (length points) do
	   (let ((dist2 (norm2 (v- (points i) point))))
	     (when (< dist2 nearest-dist2)
	       (setf nearest i
		     nearest-dist2 dist2))))
      nearest)))
;; 100000 2D points take 62 ms doesn't matter if component is small
;; integer or double, 3D points take 200 ms
#+nil
(let* ((n 100000)
	(a (make-array n :element-type 'vec
		       :initial-contents (loop repeat n collect
					      (v (random 20d0) (random 20d0)))))
	(v (v 5d0 5d0 10d0)))
   (time (nearest-naive a v)))


(deftype axis ()
  `(member ,@(loop for i below (1- +k+) collect i)))

(defun next-axis (axis)
  (declare (axis axis))
  (mod (1+ axis) +k+))

;;   (min (v (- +big+) (- +big+) (- +big+)) :type vec)
;;   (max (v +big+ +big+ +big+) :type vec)

(defstruct node 
  (point (required-argument :point) :type vec)
  (left nil :type (or null node))
  (right nil :type (or null node)))

#+nil (defstruct (interior (:include node))
  (left (required-argument :left) :type node)
  (right (required-argument :right) :type node)
  (axis (required-argument :axis) :type axis))

#+nil(defstruct (leaf (:include node))
  )

#+nil(defun in-region-p (region-min region-max point)
  (declare (vec region-min region-max point)
	   (values boolean &optional))
  (with-arrays (region-min region-max point)
    (let ((eps 1d-12)) ;; account for roundoff error in position
     (dotimes (i +k+)
       (when (or (< (point i) (- (region-min i) eps))
		 (< (+ (region-max i) eps) (point i)))
	 (return-from in-region-p nil))))
   t))
#+nil
(let ((mi (v))
      (ma (v 10.0 10.0))
      (p2 (v 10.0 (* 2.0 5.0)))
      (p (v 1.0 1.0))
      (q (v 11.0 11.0)))
  (list
   (in-region-p mi ma p)
   (in-region-p mi ma q)
   (in-region-p mi mi mi)
   (in-region-p ma ma ma)
   (in-region-p ma ma p2)))
;; should be t nil t t t


(defun successor (axis root q)
  (declare (axis axis)
	   (node root q)
	   (values node &optional))
  (let* ((kr (aref (node-point root) axis))
	 (kq (aref (node-point q) axis)))
    ;; do I need the superkey?
    (if (< kq kr)
	(node-left root)
	(node-right root))))

(defun insert (axis root p)
  (declare (axis axis)
	   ((or null node) root)
	   (node p))
  (let* ((q (if root
		root
		(return-from insert p)))
	 (son (successor axis root q)))
    (if son
	(insert axis son p)
	(setf son p))))
#+nil
(insert 0 nil (make-node :point (v .1d0 .3d0)))
#+nil
(let* ((n 10)
       (points (loop for i below n 
		  collect (v (random 1d0) (random 1d0))))
       (tree nil))
  (dolist (p points)
    (setf tree (insert 0 tree (make-node :point p))))
  tree)
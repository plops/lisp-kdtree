#.(require :alexandria)
(defpackage :kdtree (:use :cl :alexandria))
(in-package :kdtree)
(declaim (optimize (speed 2) (debug 3) (safety 3)))

(defconstant +k+ 2)
(defconstant +big+ 1e12)

(deftype vec ()
  `(simple-array float (,+k+)))

(defun v (&optional (x .0) (y .0) (z .0))
  (declare (float x y)
	   (values vec &optional))
  (make-array +k+ :element-type 'float
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
(do-vec (i (a b) c)
  (setf (c i) (+ (a i) (b i))))

#+nil
(do-vec (i () c)
  (setf (c i) (+ (a i) (b i))))
#+nil
(do-vec (i () )
  (setf (c i) (+ (a i) (b i))))

#+nil
(do-vec (i (b) )
  (setf (c i) (+ (a i) (b i))))

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
	   (float s)
	   (values vec &optional))
  (do-vec (i (a) c)
    (setf (c i) (* s (a i)))))

(defun v. (a b)
  (declare (vec a b)
	   (values float &optional))
  (let ((sum .0))
    (do-vec (i (a b))
      (incf sum (* (a i) (b i))))
    sum))
#+nil
(v. (v .1) (v 2.0))

(deftype axis ()
  `(member ,@(loop for i below (1- +k+) collect i)))

(defun next-axis (axis)
  (declare (axis axis))
  (mod (1+ axis) +k+))

(defstruct node 
  (min (v (- +big+) (- +big+) (- +big+)) :type vec)
  (max (v +big+ +big+ +big+) :type vec))

(defstruct (interior (:include node))
  (left (required-argument :left) :type node)
  (right (required-argument :right) :type node)
  (axis (required-argument :axis) :type axis))

(defstruct (leaf (:include node))
  (point (required-argument) :type vec))
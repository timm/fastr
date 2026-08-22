; vim: set ft=lisp ts=2 sw=2 et :
; lib.lisp -- batteries: macros, random, strings, csv, cli.
; (c) 2026 Tim Menzies <timm@ieee.org> MIT license
#+sbcl (declaim (sb-ext:muffle-conditions
                 warning style-warning))
(setf *read-default-float-format* 'double-float)

(defvar *the* nil) ; app settings; an app must set this

;;; ----- macros ------------------------------------------------
(defmacro ? (k) ; setf-able settings access: (? bins)
  `(,k *the*))

(defmacro fn (&body b) ; short lambda: % %1 %2; %* = all args
  `(lambda (&rest %*)
     (declare (ignorable %*))
     (symbol-macrolet
         ((% (nth 0 %*)) (%1 (nth 1 %*)) (%2 (nth 2 %*)))
       ,@b)))

(defmacro place (name &rest a) ; one spec --> get + set
  `(progn
     (defun ,name (x k &optional d)
       (declare (ignorable d))
       (typecase x ,@a))
     (defun (setf ,name) (v x k &optional d)
       (declare (ignorable d))
       (typecase x 
         ,@(loop for (ty f) in a collect `(,ty (setf ,f v)))))))

(place at ; one accessor, any container, setf-able
  (hash-table (gethash k x d))
  (vector     (aref x k))
  (list       (nth k x))
  (t          (slot-value x k)))

(set-macro-character #\$ ; $x reads as (at i 'x); i = self
  (lambda (s c) (declare (ignore c)) `(at i ',(read s t nil t)))
  t)

;;; ----- random ------------------------------------------------
(defun rand (&optional (n 1))
  (setf (? seed) (mod (* 16807 (? seed)) 2147483647))
  (* n (/ (? seed) 2147483647d0)))

(defun rint (n) (floor (rand n)))

;;; ----- strings, csv ------------------------------------------
(defun thing (s)
  (let* ((s (string-trim " " s))
         (x (let ((*read-eval* nil)) (read-from-string s nil))))
    (if (numberp x) x s)))

(defun cells (s &optional (a 0)) ; "1,2" --> (1 2)
  (let ((b (position #\, s :start a)))
    (cons (thing (subseq s a b))
          (and b (cells s (1+ b))))))

(defun csv (file fn) ; call fun on each line's cells
  (with-open-file (s file)
    (loop (funcall fn (cells (or (read-line s nil) (return)))))))

;;; ----- cli ---------------------------------------------------
(defun main (&aux (seed0 (? seed)))
  (loop for (s v) on (cdr sb-ext:*posix-argv*) do
    (let* ((w (intern (string-upcase (string-left-trim "-" s))))
           (f (intern (format nil "EG-~a" w))))
      (cond ((fboundp f)
             (setf (? seed) seed0)
             (funcall f))
            ((slot-exists-p *the* w)
             (setf (at *the* w) (thing v))
             (when (eq w 'seed)
               (setf seed0 (? seed))))))))

;;; ----- demos -------------------------------------------------
(defun eg-at () ; one accessor: hash, vector, list, struct
  (let ((h (make-hash-table :test #'equal))
        (v (vector 1 2 3)) (l (list 10 20 30)))
    (setf (at h "k") 5 (at v 1) 22 (at l 2) 33)
    (incf (at h "k"))
    (incf (at h 'fred 0) 10)
    (assert (and (= (at h "k") 6) (= (at v 1) 22)
                 (= (at l 2) 33)  (= (at h 'fred) 10)))
    (format t "hash ~a vec ~a list ~a fred ~a~%"
            (at h "k") (at v 1) (at l 2) (at h 'fred))))

(defun eg-csv () ; csv calls a function on each row
  (let ((n 0) first)
    (csv (? file) (fn (unless first (setf first %))
                      (incf n)))
    (assert (> n 100))
    (format t "~a rows; first ~a~%" n first)))

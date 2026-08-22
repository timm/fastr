; vim: set ft=lisp ts=2 sw=2 et :
; tiny.lisp -- super, smallest: cluster, contrast, classify,
; regress, multi-objective explanations.
; (c) 2026 Tim Menzies <timm@ieee.org> MIT license
(load (merge-pathnames "lib.lisp" *load-truename*))

(defstruct (settings (:conc-name nil))
  (p 2) (stop 4) (few 256) (k 5) (check 5) (bins 10)
  (cliffs 0.197) (cohen 0.35) (seed 1234) (big 1d32)
  (file "/Users/timm/gits/moot/optimize/misc/auto93.csv"))

(setf *the* (make-settings))

;;; ----- columns ----------------------------------------------
(defstruct (num (:conc-name) (:constructor num ()))
  (n 0) (mu 0d0) (m2 0d0))

(defun sym () (make-hash-table :test #'equal))

(defmethod add ((i num) v &optional (w 1)) ; welford
  (unless (equal v "?")
    (incf $n w)
    (let ((d (- v $mu)))
      (incf $mu (/ (* w d) $n))
      (incf $m2 (* w d (- v $mu)))))
  i)

(defmethod add ((i hash-table) v &optional (w 1))
  (unless (equal v "?") (incf (at i v 0) w))
  i)

(defun adds (l &optional (i (num)))
  (dolist (v l i) (add i v)))

(defun mid (i)
  (if (num-p i) $mu
      (let (mode (hi -1))
        (maphash (fn (when (> %1 hi) (setf hi %1 mode %)))
                 i)
        mode)))

(defun sd (i)
  (if (< $n 2) 0 (sqrt (/ (max 0 $m2) (1- $n)))))

(defun ent (i)
  (let ((all 0) (e 0))
    (maphash (fn (incf all %1)) i)
    (maphash (fn (when (> %1 0)
                   (decf e (* (/ %1 all)
                              (log (/ %1 all) 2)))))
             i)
    e))

(defun norm (i v)
  (if (equal v "?") v
      (let ((z (/ (- v $mu) (+ (sd i) (/ 1 (? big))))))
        (/ 1 (+ 1 (exp (* -1.7 (max -3 (min 3 z)))))))))

;;; ----- demos ------------------------------------------------
(defun eg-num () ; Num tracks n, mu, sd
  (let ((i (num)))
    (dotimes (j 100) (add i (rand)))
    (assert (and (= $n 100) (< (abs (- $mu .5)) .1)))
    (assert (< (abs (- (sd i) (/ 1 (sqrt 12)))) .05))
    (format t "n ~a mu ~,3f sd ~,3f~%" $n $mu (sd i))))

(defun eg-sym () ; Sym counts; mid is mode; ent is spread
  (let ((i (adds '("a" "a" "b" "b" "b" "c") (sym))))
    (assert (and (= (at i "b") 3) (equal (mid i) "b")))
    (format t "b ~a mid ~a ent ~,2f~%"
            (at i "b") (mid i) (ent i))))

(main)

; unsuper: unsupervised recursive bi-clustering (no goals).
; Fastmap poles halve the data; contrast names the split.
; (c) 2026 Tim Menzies <timm@ieee.org> MIT license
(defvar *loading* nil)
(let ((*loading* t))
  (load (merge-pathnames "contrast.lisp" *load-truename*)))

(defstruct node rows n cut yes no ys) ; one tree node

(defun distx (tbl r1 r2) ; distance between non-goal values
  (let ((d 0) (n 1d-32))
    (dolist (at (cols-x (tbl-cols tbl)))
      (incf n)
      (incf d (expt (distx1 (col-of tbl at)
                            (nth at r1) (nth at r2))
                    (? :p))))
    (expt (/ d n) (/ 1 (? :p)))))

(defun distx1 (col a b) ; helper for one column
  (cond ((and (equal a "?") (equal b "?")) 1)
        ((not (num-p col)) (if (equal a b) 0 1))
        (t (let ((a (norm col a)) (b (norm col b)))
             (when (equal a "?")
               (setf a (if (< b .5) 1 0)))
             (when (equal b "?")
               (setf b (if (< a .5) 1 0)))
             (abs (- a b))))))

(defun proj (tbl lohi row) ; project row onto line lo -> hi
  (let ((a (distx tbl row (at lohi 'lo)))
        (b (distx tbl row (at lohi 'hi)))
        (c (at lohi 'c)))
    (/ (+ (* a a) (* c c) (- (* b b))) (+ (* 2 c) 1d-32))))

(defun fastmap (tbl rows) ; sort rows on line between 2 poles
  (let* ((some (many rows (min (? :few) (length rows))))
         (w  (any some))
         (m1 (most (fn (distx tbl w %)))))
    (map nil m1 some)
    (let* ((a (funcall m1))
           (m2 (most (fn (distx tbl a %)))))
      (map nil m2 some)
      (let* ((b (funcall m2))
             (lohi (list (cons 'lo a) (cons 'hi b)
                         (cons 'c (distx tbl a b)))))
        (values lohi
                (keysort rows (fn (proj tbl lohi %))))))))

(defun cluster (tbl &optional rows) ; contrast names the splits
  (labels
    ((go1 (rows)
       (let ((tr (make-node :rows rows :n (length rows))))
         (when (> (length rows) (? :stop))
           (let* ((tmp (nth-value 1 (fastmap tbl rows)))
                  (half (floor (length tmp) 2))
                  (cut (contrasts tbl (subseq tmp 0 half)
                                      (subseq tmp half))))
             (when cut
               (let ((yes (remove-if-not
                            (fn (selects cut %)) rows))
                     (no  (remove-if
                            (fn (selects cut %)) rows)))
                 (when (< 0 (length yes) (length rows))
                   (setf (node-cut tr) cut
                         (node-yes tr) (go1 yes)
                         (node-no tr)  (go1 no)))))))
         tr)))
    (go1 (or rows (tbl-rows tbl)))))

(defun show (tr &optional pre txt more) ; print tree with n
  (format t "~a~an=~a~a~%" (or pre "") (or txt "")
          (node-n tr)
          (if more (funcall more tr) ""))
  (when (node-cut tr)
    (let ((sub (if (null pre) "" (concatenate 'string
                                   pre "|  ")))
          (s (at (node-cut tr) 'txt)))
      (show (node-yes tr) sub
            (format nil "~a; " s) more)
      (show (node-no tr) sub
            (format nil "!~a; " s) more))))

(defun leaves (tr) ; all terminal nodes of a tree
  (if (node-cut tr)
      (append (leaves (node-yes tr)) (leaves (node-no tr)))
      (list tr)))

;;; ------------------------------------------------------------
(defun test-distx () ; symmetric, zero on self, in 0..1
  (let* ((tb (tbl+ (? :file)))
         (rows (tbl-rows tb))
         (r1 (aref rows 0))
         (r2 (aref rows (1- (length rows)))))
    (assert (zerop (distx tb r1 r1)))
    (assert (= (distx tb r1 r2) (distx tb r2 r1)))
    (loop for r across rows do
      (assert (<= 0 (distx tb r1 r) 1)))
    (format t "d(r1,r1) ~,2f | d(r1,r2) ~,2f~%"
            (distx tb r1 r1) (distx tb r1 r2))
    (format t "r1 ~a~%r2 ~a~%" r1 r2)))

(defun test-fastmap () ; rows sorted along a line of 2 poles
  (let ((tb (tbl+ (? :file))))
    (multiple-value-bind (lohi rows) (fastmap tb (tbl-rows tb))
      (assert (> (at lohi 'c) 0))
      (let ((last -1d32))
        (loop for r across rows
              for p = (proj tb lohi r) do
          (assert (>= p last)) (setf last p)))
      (format t "pole1 ~a~%pole2 ~a~%separation ~,2f~%"
              (at lohi 'lo) (at lohi 'hi) (at lohi 'c)))))

(defun test-cluster () ; print tree; leaf counts sum to n
  (let* ((tb (tbl+ (? :file)))
         (tr (cluster tb)))
    (assert (= (node-n tr) (length (tbl-rows tb))))
    (assert (= (node-n tr)
               (reduce #'+ (mapcar #'node-n (leaves tr)))))
    (show tr)
    (format t "~a leaves~%" (length (leaves tr)))))

(unless *loading* (main *load-truename*))

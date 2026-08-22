; contrast: ranges that most separate two groups of rows.
; X ranges come from the parent table; one XY per column counts
; which group lands in each range. Contrast is just counting.
; (c) 2026 Tim Menzies <timm@ieee.org> MIT license
(defvar *loading* nil)
(let ((*loading* t))
  (load (merge-pathnames "tbl.lisp" *load-truename*)))

(defstruct span lo hi first last score at txt)

(defun span+ (lo hi b r first last) ; the scorer: b^2/(b+r)
  (make-span :lo lo :hi hi :first first :last last
             :score (/ (* b b) (+ b r 1d-32))))

(defun discretize (d0 d1 ordered out) ; two count hashes
  (let* ((n0 (+ 1d-32 (loop for v being the hash-values of d0
                            sum v)))
         (n1 (+ 1d-32 (loop for v being the hash-values of d1
                            sum v)))
         (ks (sort (union
                     (loop for k being the hash-keys of d0
                           collect k)
                     (loop for k being the hash-keys of d1
                           collect k) :test #'equal)
                   (fn (if (numberp %) (< % %1) t))))
         (m (length ks))
         (b+ (list 0)) (r+ (list 0)))
    (dolist (k ks)
      (push (+ (first b+) (/ (gethash k d0 0) n0)) b+)
      (push (+ (first r+) (/ (gethash k d1 0) n1)) r+))
    (let ((b+ (coerce (nreverse b+) 'vector))
          (r+ (coerce (nreverse r+) 'vector))
          (ks (coerce ks 'vector)))
      (dotimes (i m)
        (loop for j from (1+ i)
              to (if ordered m (1+ i)) do
          (unless (and (= i 0) (= j m)) ; whole col useless
            (funcall out
              (span+ (aref ks i) (aref ks (1- j))
                     (- (aref b+ j) (aref b+ i))
                     (- (aref r+ j) (aref r+ i))
                     (= i 0) (= j m)))))))))

(defun contrasts (tbl rows0 rows1) ; best range splitting groups
  (let ((cs (ranges tbl))
        (all (most #'span-score)))
    (dolist (at (cols-x (tbl-cols tbl)) (funcall all))
      (let ((cuts (gethash at cs))
            (xy (make-xy :cuts (gethash at cs)))
            (one (most #'span-score)))
        (loop for (y rows) in (list (list 0 rows0)
                                    (list 1 rows1)) do
          (map nil (fn (add xy (list (nth at %) y))) rows))
        (let ((d0 (sym+)) (d1 (sym+)))
          (loop for b being the hash-keys of (xy-seen xy)
                using (hash-value d) do
            (setf (gethash b d0) (gethash 0 d 0)
                  (gethash b d1) (gethash 1 d 0)))
          (discretize d0 d1 cuts one)
          (discretize d1 d0 cuts one))
        (let ((z (funcall one)))
          (when z
            (setf (span-at z) at)
            (let ((name (name-of tbl at)))
              (if (null cuts)
                  (setf (span-txt z)
                        (format nil "~a == ~a" name
                                (span-lo z)))
                  (progn
                    (setf (span-lo z)
                          (if (span-first z) -1d32
                              (nth (span-lo z) cuts))
                          (span-hi z)
                          (if (span-last z) 1d32
                              (nth (1+ (span-hi z)) cuts)))
                    (setf (span-txt z)
                      (cond ((= (span-lo z) -1d32)
                             (format nil "~a <= ~a" name
                                     (g3 (span-hi z))))
                            ((= (span-hi z) 1d32)
                             (format nil "~a >= ~a" name
                                     (g3 (span-lo z))))
                            (t (format nil "~a in ~a..~a"
                                       name (g3 (span-lo z))
                                       (g3 (span-hi z)))))))))
            (funcall all z)))))))

(defun selects (cut row) ; does row fall inside the range?
  (let ((v (nth (span-at cut) row)))
    (and (not (equal v "?"))
         (if (numberp (span-lo cut))
             (and (numberp v)
                  (<= (span-lo cut) v (span-hi cut)))
             (equal v (span-lo cut))))))

;;; ------------------------------------------------------------
(defun test-discretize () ; separated count dicts; span covers d0
  (let ((d0 (sym+)) (d1 (sym+)) (out (most #'span-score)))
    (loop for (k v) in '((0 30) (1 50) (2 20)) do
      (setf (gethash k d0) v))
    (loop for (k v) in '((6 40) (7 40) (8 20)) do
      (setf (gethash k d1) v))
    (discretize d0 d1 t out)
    (let ((z (funcall out)))
      (assert (and (> (span-score z) .9)
                   (= (span-lo z) 0) (= (span-hi z) 2)))
      (format t "d0 bins 0-2, d1 bins 6-8: span ~a..~a ~,2f~%"
              (span-lo z) (span-hi z) (span-score z)))))

(defun test-contrasts () ; name fastmap's 50:50 top split
  (let ((*loading* t))
    (load (merge-pathnames "unsuper.lisp" *load-truename*)))
  (let* ((tb (tbl+ (? :file)))
         (rows (nth-value 1 (fastmap tb (tbl-rows tb))))
         (half (floor (length rows) 2))
         (r0 (subseq rows 0 half)) (r1 (subseq rows half))
         (z (contrasts tb r0 r1))
         (p0 (/ (count-if (fn (selects z %)) r0)
                (length r0)))
         (p1 (/ (count-if (fn (selects z %)) r1)
                (length r1))))
    (format t "best: ~a; selects ~,0f% vs ~,0f%~%"
            (span-txt z) (* 100 p0) (* 100 p1))
    (assert (> (abs (- p0 p1)) .2))))

(defun test-halve () ; one line: best split of two halves
  (let ((*loading* t))
    (load (merge-pathnames "unsuper.lisp" *load-truename*)))
  (let* ((tb (tbl+ (? :file)))
         (rows (nth-value 1 (fastmap tb (tbl-rows tb))))
         (half (floor (length rows) 2))
         (z (contrasts tb (subseq rows 0 half)
                          (subseq rows half))))
    (format t "~45a ~,2f ~a~%"
            (file-namestring (? :file))
            (if z (span-score z) 0)
            (if z (span-txt z) "no split"))))

(unless *loading* (main *load-truename*))

; super: supervised recursive splitting for classify, regress.
; Each split minimizes the spread of the goals in two halves.
; (c) 2026 Tim Menzies <timm@ieee.org> MIT license
(defvar *loading* nil)
(let ((*loading* t))
  (load (merge-pathnames "unsuper.lisp" *load-truename*)))

(defun disty (tbl row) ; distance of goals to best corner
  (let ((d 0) (n 1d-32))
    (loop for (at . w) in (cols-y (tbl-cols tbl))
          for v = (norm (col-of tbl at) (nth at row))
          unless (equal v "?") do
            (incf n) (incf d (expt (abs (- v w)) (? :p))))
    (expt (/ d n) (/ 1 (? :p)))))

(defun yinfo (tbl) ; goal fun, summary maker, spread fun
  (if (cols-y (tbl-cols tbl))
      (values (fn (disty tbl %)) #'num+ #'sd)
      (let ((kl (cols-klass (tbl-cols tbl))))
        (values (fn (nth kl %)) #'sym+ #'ent))))

(defun cut+ (tbl at cut score) ; a supervised cut record
  (list (cons 'at at) (cons 'cut cut) (cons 'score score)
        (cons 'txt (format nil "~a <= ~a"
                           (name-of tbl at) (g3 cut)))))

(defun split (tbl rows) ; col+cut minimizing goal spread
  (multiple-value-bind (y new div) (yinfo tbl)
    (let (out)
      (dolist (at (cols-x (tbl-cols tbl)) out)
        (when (num-p (col-of tbl at))
          (setf out (split1 tbl rows at y new div out)))))))

(defun split1 (tbl rows at y new div out) ; try one column
  (let ((have (coerce (keysort
                        (remove-if
                          (fn (equal (nth at %) "?")) rows)
                        (fn (nth at %)))
                      'vector)))
    (when (>= (length have) 4)
      (let* ((ys (map 'vector y have))
             (left (funcall new))
             (right (adds ys (funcall new)))
             (m (length have)))
        (dotimes (i (1- m))
          (add left (aref ys i))
          (add right (aref ys i) -1)
          (let ((v0 (nth at (aref have i)))
                (v1 (nth at (aref have (1+ i)))))
            (unless (= v0 v1)
              (let ((s (/ (+ (* (1+ i) (funcall div left))
                             (* (- m i 1)
                                (funcall div right)))
                          m)))
                (when (or (null out) (< s (at out 'score)))
                  (setf out (cut+ tbl at
                                  (/ (+ v0 v1) 2) s)))))))))
    out))

(defun tree (tbl &optional rows) ; recursive supervised splits
  (multiple-value-bind (y new) (yinfo tbl)
    (labels
      ((below (cut row)
         (let ((v (nth (at cut 'at) row)))
           (and (not (equal v "?")) (<= v (at cut 'cut)))))
       (go1 (rows)
         (let ((tr (make-node
                     :rows rows :n (length rows)
                     :ys (adds (map 'vector y rows)
                               (funcall new)))))
           (when (> (length rows) (? :stop))
             (let ((cut (split tbl rows)))
               (when cut
                 (let ((yes (remove-if-not
                              (fn (below cut %)) rows))
                       (no  (remove-if
                              (fn (below cut %)) rows)))
                   (when (< 0 (length yes) (length rows))
                     (setf (node-cut tr) cut
                           (node-yes tr) (go1 yes)
                           (node-no tr)  (go1 no)))))))
           tr)))
      (go1 (or rows (tbl-rows tbl))))))

(defun leaf (tr row) ; walk row down to its leaf
  (loop while (node-cut tr)
        for cut = (node-cut tr)
        for v = (nth (at cut 'at) row)
        do (setf tr (if (and (not (equal v "?"))
                             (<= v (at cut 'cut)))
                        (node-yes tr) (node-no tr))))
  tr)

(defun predict (tbl tr row &optional (k (? :k))) ; knn in leaf
  (let* ((rows (keysort (node-rows (leaf tr row))
                        (fn (distx tbl row %))))
         (kl (cols-klass (tbl-cols tbl)))
         (got (col+ (name-of tbl kl))))
    (loop for i below (min k (length rows)) do
      (add got (nth kl (aref rows i))))
    (mid got)))

;;; ------------------------------------------------------------
(defun test-disty () ; sort rows by goals; show best, worst
  (let* ((tb (tbl+ (? :file)))
         (lst (keysort (tbl-rows tb) (fn (disty tb %))))
         (m (length lst)))
    (assert (< (disty tb (aref lst 0))
               (disty tb (aref lst (1- m)))))
    (format t "~a~%" (cols-names (tbl-cols tb)))
    (dolist (j (list 0 1 2 (- m 3) (- m 2) (- m 1)))
      (format t "~a d ~,2f~%" (aref lst j)
              (disty tb (aref lst j))))))

(defun test-split () ; show split minimizing goal spread
  (let* ((tb (tbl+ (? :file)))
         (cut (split tb (tbl-rows tb))))
    (format t "~a spread ~,3f~%"
            (at cut 'txt) (at cut 'score))))

(defun test-tree () ; print tree with n and mean goal dist
  (let* ((tb (tbl+ (? :file)))
         (tr (tree tb)))
    (assert (= (node-n tr)
               (reduce #'+ (mapcar #'node-n (leaves tr)))))
    (show tr nil nil
          (fn (format nil " mu=~,2f" (mid (node-ys %)))))
    (format t "~a leaves~%" (length (leaves tr)))))

(defun test-classify () ; resubstitution accuracy, diabetes
  (let* ((tb (tbl+
               "/Users/timm/gits/moot/classify/diabetes.csv"))
         (tr (tree tb))
         (kl (cols-klass (tbl-cols tb)))
         (acc (/ (count-if
                   (fn (equal (predict tb tr %) (nth kl %)))
                   (tbl-rows tb))
                 (length (tbl-rows tb)))))
    (format t "resub acc ~,2f~%" (float acc))
    (assert (> acc .7))))

(unless *loading* (main *load-truename*))

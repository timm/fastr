; tbl: Num, Sym, Cols, Tbl, XY: incremental column summaries.
; (c) 2026 Tim Menzies <timm@ieee.org> MIT license
(defvar *loading* nil)
(let ((*loading* t))
  (load (merge-pathnames "start.lisp" *load-truename*)))

(defstruct (num (:constructor num+ ())) ; n mu m2 of numbers
  (n 0) (mu 0d0) (m2 0d0))

(defun sym+ () ; counts of symbols
  (make-hash-table :test #'equal))

(defstruct cols names all x y klass) ; columns grouped x,y

(defstruct tbl ; rows plus their column summaries
  (rows (make-array 0 :adjustable t :fill-pointer t)) cols)

(defstruct xy cuts (y #'sym+) (seen (make-hash-table)))

(defun col+ (name) ; Nums start with an upper case letter
  (if (upper-case-p (char name 0)) (num+) (sym+)))

(defun cols+ (names) ; names --> columns grouped into x,y
  (let ((i (make-cols :names names)))
    (loop for name in names for at from 0
          for z = (char name (1- (length name))) do
      (push (col+ name) (cols-all i))
      (cond ((char= z #\X))
            ((char= z #\!) (setf (cols-klass i) at))
            ((char= z #\+) (push (cons at 1) (cols-y i)))
            ((char= z #\-) (push (cons at 0) (cols-y i)))
            (t             (push at (cols-x i)))))
    (setf (cols-all i) (nreverse (cols-all i))
          (cols-x i)   (nreverse (cols-x i))
          (cols-y i)   (nreverse (cols-y i)))
    i))

(defun tbl+ (file &optional i) ; csv file --> table
  (csv file (fn (if i (add i %)
                    (setf i (make-tbl :cols (cols+ %))))))
  i)

(defun clone (tbl &optional rows) ; empty copy, same columns
  (let ((i (make-tbl :cols (cols+ (cols-names (tbl-cols tbl))))))
    (map nil (fn (add i %)) (or rows #()))
    i))

(defun col-of (tbl at) (nth at (cols-all (tbl-cols tbl))))
(defun name-of (tbl at) (nth at (cols-names (tbl-cols tbl))))

(defun add (i v &optional (w 1)) ; add value (or row) v
  (unless (equal v "?")
    (typecase i
      (num        (welford i v w))
      (hash-table (incf (gethash v i 0) w))
      (cols       (mapc (fn (add % %1 w)) (cols-all i) v))
      (tbl        (vector-push-extend v (tbl-rows i))
                  (add (tbl-cols i) v w))
      (xy         (destructuring-bind (x y) v
                    (unless (or (equal x "?") (equal y "?"))
                      (let* ((b (if (xy-cuts i)
                                    (bin (xy-cuts i) x) x))
                             (seen (xy-seen i)))
                        (setf (gethash b seen)
                              (add (or (gethash b seen)
                                       (funcall (xy-y i)))
                                   y w))))))))
  i)

(defun adds (seq &optional (i (num+))) ; add all of seq into i
  (map nil (fn (add i %)) seq)
  i)

(defun welford (i v w) ; incrementally update a Num
  (with-slots (n mu m2) i
    (incf n w)
    (if (< n 1) (setf n 0 mu 0d0 m2 0d0)
        (let ((d (- v mu)))
          (incf mu (/ (* w d) n))
          (incf m2 (* w d (- v mu)))))))

(defun mid (col) ; middle of a distribution
  (if (num-p col) (num-mu col)
      (let (mode (hi -1))
        (loop for k being the hash-keys of col
              using (hash-value n)
              when (> n hi) do (setf hi n mode k))
        mode)))

(defun mids (tbl) ; middles of all columns
  (mapcar #'mid (cols-all (tbl-cols tbl))))

(defun sd (num) ; diversity of a Num
  (with-slots (n m2) num
    (if (<= n 1) 0 (sqrt (/ m2 (- n 1))))))

(defun ent (d) ; diversity of a Sym
  (let ((all 0) (e 0))
    (loop for n being the hash-values of d do (incf all n))
    (loop for n being the hash-values of d
          when (> n 0) do
            (decf e (* (/ n all) (log (/ n all) 2))))
    e))

(defun norm (col x) ; x --> 0..1 via logistic cdf approx
  (if (equal x "?") x
      (let ((z (/ (- x (num-mu col)) (+ (sd col) 1d-32))))
        (/ 1 (+ 1 (exp (* -1.7 (max -3 (min 3 z)))))))))

(defun ranges (tbl) ; per numeric x col: bins cuts, +/-3sd
  (let ((cs (make-hash-table)))
    (dolist (at (cols-x (tbl-cols tbl)) cs)
      (let ((col (col-of tbl at)))
        (when (num-p col)
          (let* ((lo (- (num-mu col) (* 3 (sd col))))
                 (hi (+ (num-mu col) (* 3 (sd col)))))
            (setf (gethash at cs)
                  (loop for k to (? :bins) collect
                        (+ lo (/ (* k (- hi lo))
                                 (? :bins)))))))))))

(defun bin (cuts x) ; index of the cut interval holding x
  (let ((lo (first cuts)) (hi (car (last cuts))))
    (max 0 (min (- (length cuts) 2)
                (floor (/ (* (- x lo) (1- (length cuts)))
                          (+ (- hi lo) 1d-32)))))))

;;; ------------------------------------------------------------
(defun test-num () ; Num tracks n,mu,sd; subtract restores
  (let ((i (num+)))
    (loop repeat 100 do (add i (rand)))
    (let ((n (num-n i)) (mu (num-mu i)))
      (assert (and (= n 100) (< (abs (- mu .5)) .1)))
      (assert (< (abs (- (sd i) (/ 1 (sqrt 12)))) .05))
      (add (add i .5) .5 -1)
      (assert (and (< (abs (- (num-mu i) mu)) 1d-9)
                   (= (num-n i) n)))
      (format t "100 rands: n ~a mu ~,3f sd ~,3f~%"
              n mu (sd i)))))

(defun test-sym () ; Sym counts symbols; mid is the mode
  (let ((i (sym+)))
    (loop for c across "aabbbc" do (add i (string c)))
    (assert (and (= (gethash "b" i) 3) (equal (mid i) "b")))
    (format t "aabbbc: b=~a mid=~a~%" (gethash "b" i) (mid i))))

(defun test-tbl () ; load file; find rows and x,y columns
  (let ((tb (tbl+ (? :file))))
    (assert (> (length (tbl-rows tb)) 100))
    (format t "~a rows; x: ~a y: ~a~%"
            (length (tbl-rows tb))
            (cols-x (tbl-cols tb)) (cols-y (tbl-cols tb)))))

(defun test-mid () ; show middle of every column
  (let ((tb (tbl+ (? :file))))
    (format t "~{~,5@a ~}~%"
            (mapcar (fn (if (floatp %) (format nil "~,1f" %) %))
                    (mids tb)))))

(defun test-norm () ; norm maps to 0..1, monotonically
  (let ((i (num+)) (last 0))
    (loop repeat 100 do (add i (rand 0 1000)))
    (loop for x from 0 to 999 by 99
          for v = (norm i x) do
            (assert (and (<= last v) (<= 0 v 1)))
            (setf last v))
    (loop for x in '(0 250 500 750 999) do
      (format t "~a->~,2f " x (norm i x)))
    (terpri)))

(defun test-cdf () ; norm(hi)-norm(lo) ~ mass inside lo..hi
  (let ((i (num+)) xs)
    (loop repeat 1000 do
      (let ((x (- (loop repeat 12 sum (rand)) 6)))
        (push x xs) (add i x)))
    (loop for (lo hi) in '((-1 1) (0 1) (-2 2)) do
      (let ((got  (- (norm i hi) (norm i lo)))
            (want (/ (count-if (fn (<= lo % hi)) xs)
                     (length xs))))
        (format t "mass ~2a..~a: cdf ~,2f truth ~,2f~%"
                lo hi got want)
        (assert (< (abs (- got want)) .05))))))

(defun test-xy () ; per x-range: y counts (Sym) or mu,sd (Num)
  (let* ((tb (tbl+ (? :file)))
         (cs (gethash 0 (ranges tb)))
         (s (position "origin" (cols-names (tbl-cols tb))
                      :test #'equal))
         (n (position "Lbs-" (cols-names (tbl-cols tb))
                      :test #'equal))
         (xy1 (make-xy :cuts cs))
         (xy2 (make-xy :cuts cs :y #'num+)))
    (loop for r across (tbl-rows tb) do
      (add xy1 (list (nth 0 r) (nth s r)))
      (add xy2 (list (nth 0 r) (nth n r))))
    (loop for b in (sort (loop for k being the hash-keys
                               of (xy-seen xy1) collect k) #'<)
          for d = (gethash b (xy-seen xy1))
          for y = (gethash b (xy-seen xy2)) do
      (format t "bin ~2a origins ~a Lbs mu ~5,0f sd ~4,0f~%"
              b (loop for k being the hash-keys of d
                      using (hash-value v) collect
                      (list k v))
              (num-mu y) (sd y)))))

(unless *loading* (main *load-truename*))

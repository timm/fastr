; start: config, atoms, csv rows, and a tiny test-runner CLI.
; (c) 2026 Tim Menzies <timm@ieee.org> MIT license
#+sbcl (declaim (sb-ext:muffle-conditions
                 warning style-warning))

(defvar *help* "
start: config, atoms, csv rows, and a tiny test-runner CLI.
(c) 2026 Tim Menzies <timm@ieee.org> MIT license

Options:

  -h             show help
  --p=2          minkowski coefficient
  --stop=4       stopping rule for recursive tree generation
  --few=256      sub-sample size for pole finding
  --k=5          nearest neighbors used in a leaf
  --check=5      optimization: how many top picks to evaluate
  --bins=10      number of bins for discretization
  --cliffs=.197  cliffs delta: max effect size for same
  --cohen=.35    cohen d: max mean separation for same
  --seed=1234    random number generation
  --file=/Users/timm/gits/moot/optimize/misc/auto93.csv")

;;; ----- macros ------------------------------------------------
(defmacro fn (&body b) ; short lambda: % %1 %2 = args 0 1 2
  (let ((a (gensym)))
    `(lambda (&rest ,a)
       (declare (ignorable ,a))
       (symbol-macrolet
           ((% (nth 0 ,a)) (%1 (nth 1 ,a)) (%2 (nth 2 ,a)))
         ,@b))))

(defmacro ? (k) ; settings accessor: (? :bins)
  `(gethash ,k *the*))

(defun at (x k &optional default) ; get k from any container
  (typecase x
    (hash-table (gethash k x default))
    (list (if (consp (first x))
              (let ((p (assoc k x :test #'equal)))
                (if p (cdr p) default))       ; alist
              (getf x k default)))            ; plist
    (t (if (slot-exists-p x k)
           (slot-value x k) default))))

(defun (setf at) (v x k &optional default) ; (setf (at x k) v)
  (declare (ignorable default))
  (typecase x
    (hash-table (setf (gethash k x) v))
    (list (let ((p (assoc k x :test #'equal)))
            (if p (setf (cdr p) v)
                (error "at: can't grow list ~a" x))))
    (t (setf (slot-value x k) v))))

;;; ----- atoms, strings ----------------------------------------
(defun g3 (x) ; number --> short string, 3 significant digits
  (if (integerp x) (format nil "~a" x)
      (string-right-trim " ." (format nil "~,3g" x))))

(defun thing (s) ; string --> number or trimmed string
  (let* ((s (string-trim " " s))
         (n (ignore-errors (read-from-string s))))
    (if (numberp n) n s)))

(defun split (s &optional (sep #\,)) ; string --> substrings
  (loop for i = 0 then (1+ j)
        for j = (position sep s :start i)
        collect (subseq s i j) while j))

(defvar *the* (make-hash-table))
(loop for line in (split *help* #\Newline)
      for a = (search "--" line)
      for b = (and a (position #\= line))
      for c = (and b (position #\Space line :start b))
      when b do (setf (gethash (intern (string-upcase
                       (subseq line (+ a 2) b)) :keyword) *the*)
                      (thing (subseq line (1+ b) c))))

;;; ----- random ------------------------------------------------
(defvar *seed* 1234)

(defun rand (&optional (lo 0) (hi 1)) ; pseudo-random lo..hi
  (setf *seed* (mod (* 16807 *seed*) 2147483647))
  (+ lo (* (- hi lo) (/ *seed* 2147483647d0))))

(defun rint (lo hi) ; pseudo-random integer lo..hi
  (floor (+ 0.5 (rand lo hi))))

(defun any (v) ; pick one item of vector v at random
  (aref v (rint 0 (1- (length v)))))

(defun many (v n) ; pick n items of vector v at random
  (coerce (loop repeat n collect (any v)) 'vector))

(defun shuffle (v) ; randomize a vector, in place
  (loop for i from (1- (length v)) downto 1
        for j = (rint 0 i)
        do (rotatef (aref v i) (aref v j)))
  v)

(defun keysort (v fun) ; sorted copy of vector v, order fun
  (sort (copy-seq v) #'< :key fun))

(defun most (fun) ; carried max closure: (f x) adds, (f) reads
  (let (b)
    (fn (when (and % (or (null b) (> (funcall fun %)
                                     (funcall fun b))))
          (setf b %))
        b)))

;;; ----- csv ---------------------------------------------------
(defun csv (file fun) ; call fun on each row of atoms
  (with-open-file (s file)
    (loop for line = (read-line s nil) while line
          do (funcall fun (mapcar #'thing (split line))))))

;;; ----- test runner -------------------------------------------
(defvar *loading* nil) ; true while loading a parent module

(defun loads (file here) ; load parent module, marked nested
  (let ((*loading* t))
    (load (merge-pathnames file here))))

(defun run (f) ; reseed, call f, catch crashes
  (setf *seed* (? :seed))
  (handler-case (funcall f)
    (error (e) (format t "~a~%" e))))

(defun demos (file) ; this file's test- funs, in source order
  (with-open-file (s file)
    (loop for line = (read-line s nil) while line
          for a = (search "(defun test-" line)
          for b = (and (eql a 0)
                       (position #\Space line :start 7))
          for c = (search "; " line)
          when b collect
            (list (subseq line 7 b)
                  (if c (subseq line (+ c 2)) "")))))

(defun main (file) ; single dash = demo, double dash = setting
  (let ((args (cdr sb-ext:*posix-argv*)))
    (loop while args do
      (let* ((s (pop args)))
        (cond
          ((equal s "-h")
           (format t "~a~%~%Demos:~%~%" *help*)
           (loop for (k doc) in (demos file) do
             (format t "  -~13a ~a~%" (subseq k 5) doc)))
          ((equal s "-all")
           (loop for (k doc) in (demos file) do
             (format t "~%# ~a~%" k)
             (run (symbol-function
                   (intern (string-upcase k))))))
          ((and (> (length s) 2) (string= s "--" :end1 2))
           (let* ((b (position #\= s))
                  (k (intern (string-upcase
                              (subseq s 2 b)) :keyword))
                  (v (if b (subseq s (1+ b)) (pop args))))
             (when (nth-value 1 (gethash k *the*))
               (setf (gethash k *the*) (thing v)))))
          ((char= (char s 0) #\-)
           (let ((f (ignore-errors (symbol-function
                     (intern (string-upcase
                       (format nil "test-~a"
                         (subseq s 1))))))))
             (when f (run f)))))))))

;;; ------------------------------------------------------------
(defun test-the () ; show current settings
  (loop for k being the hash-keys of *the* do
    (format t "~a=~a " k (? k)))
  (terpri))

(defun test-thing () ; strings coerce to numbers or strings
  (assert (eql (thing "2") 2))
  (assert (= (thing "2.1") 2.1))
  (assert (equal (thing " a ") "a"))
  (format t "'2' -> ~a | '2.1' -> ~a | ' a ' -> '~a'~%"
          (thing "2") (thing "2.1") (thing " a ")))

(defun test-at () ; one accessor for hashes, alists, plists
  (let ((h (make-hash-table)) (a '((x . 1))) (p '(:x 1)))
    (setf (at h :k) 5)
    (incf (at h :k 0))
    (assert (= (at h :k) 6))
    (assert (= (at a 'x) 1))
    (assert (= (at p :x) 1))
    (assert (eq (at h :zz :none) :none))
    (format t "hash ~a alist ~a plist ~a missing ~a~%"
            (at h :k) (at a 'x) (at p :x) (at h :zz :none))))

(defun test-csv () ; csv reader finds many rows in file
  (let ((n 0))
    (csv (? :file) (fn (declare (ignorable %)) (incf n)))
    (assert (> n 100))
    (format t "~a rows~%" n)))

(unless *loading* (main *load-truename*))

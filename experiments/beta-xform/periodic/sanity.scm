;
; Collection of sanity-check routines.
; Not written in c, so that I can double-check insane mistakes.
;

(define phi (* 0.5 (+ 1 (sqrt 5))))

(define beta phi)

; iterate once using the beta shift
(define (tmid x) (* beta (if (< 0.5 x) (- x 0.5) x)))

; Iterate n times
(define (tmn x n) (if (= n 0) x (tmn (tmid x) (- n 1))))

; Return the nth bit of midpoint iteration
(define (bn n) (if (< 0.5 (tmn (* 0.5 beta) n)) 1 0))

; Return the bit sequence of iterating the midpoint n times.
(define (bitseq n) (map bn (iota n)))

(define beta 1.839286755214162) ; pn=3
(define beta 1.9275619754829253) ; pn=7
(define beta 1.965948236645486) ; pn = 15

; Sloppy midpoint subdivision root finder
(define (mroot f lo hi)
	(define eps 1e-15)
	(define mid (* 0.5 (+ lo hi)))
	(define fhi (f hi))
	(define fmid (f mid))
	(if (> eps (- hi lo)) mid
		(if (> 0 (* fhi fmid)) (mroot f mid hi) (mroot f lo mid))))

(define (root f) (mroot f 1 2))

(define (oo x) (+ (* x x) (- x) -1))
(define (poz x) (+ (* x x x) (- x) 1)) 


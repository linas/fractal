
; Basic fibonacci
(define (fib n) (if (>= 1 n) 1 (+ (fib (- n 1)) (fib (- n 2)))))

; Bitseq nacci
; Does the two order-three sequences.
(define (fibit3 n b1)
	(if (>= 1 n) 1
		(+ (fibit3 (- n 1) b1) (* b1 (fibit3 (- n 2) b1)) (fibit3 (- n 3) b1))))

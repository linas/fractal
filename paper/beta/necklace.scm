

(define (neck1 a) a)
(define (neck2 a) (/ (- (* a a) a) 2))
(define (neck3 a) (/ (- (* a a a) a) 3))
(define (neck4 a) (/ (- (* a a a a ) (* a a)) 4))
(define (neck5 a) (/ (- (* a a a a a) a) 5))
(define (neck6 a) (/ (+ (- (* a a a a a a) (* a a a) (* a a)) a) 6))
(define (neck7 a) (/ (- (* a a a a a a a) a) 7))

(define (neck a)
	(format #t "1, ~D, ~D, ~D, ~D, ~D, ~D, ~D\n"
		(neck1 a)
		(neck2 a)
		(neck3 a)
		(neck4 a)
		(neck5 a)
		(neck6 a)
		(neck7 a)))



 


;; (load-module  'legendre.scm)
;;
;; compute driebes eqn (A.33) integral of rescaled legendre's
;; it defines (I m n)
;; this works for about n, m < 18 and after that,
;; start loosing accuracy 


;; define factorial function

(define (factorial n)
   (define (iter product counter max-count)
      (if (> counter max-count)
         product
         (iter (* counter product)
                    (+ counter 1)
                    max-count )))
     (iter 1 1 n))

;; x raised to power y
(define (pow x y)
    (exp (* y (log x))))

;; x raised to integer power n
(define (ipow x n)
   (define (iter product x counter max-count)
      (if (>= counter max-count)
         product
         (iter (* x product) x
                    (+ counter 1)
                    max-count )))
    (iter 1 x 0 n))
      

;; #t if n is odd
(define (odd n)  
    (= n (* 2 (floor (/ n 2)))) )

;; -1 if n odd, +1 if n even
(define (sgn n)  
    (if (odd n) 1 -1) )


(define (I m n)
   (define (term m n l)
           (/ (* (sgn l) (factorial (+ m n l)) )
            ( * (ipow 2 l)
                (factorial (- (- n m) l)) 
                (factorial (+ (* 2 m) l 1)) 
                (factorial l) ) 
            )
           )
   (define (sum m n)
      (define (iter sm m n l max-count)
         (if (> l max-count)
            sm
            (iter (+ (term m n l) sm)  
                       m n 
                       (+ l 1)
                       max-count )))
        (iter 0 m n 0 (- n m)))
   ( / (* (sgn (+ m n)) (sqrt (* (+ (* 2 m) 1) (+ (* 2 n) 1))) (sum m n))
       (ipow 2 m)))

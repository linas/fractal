
;; legendre.scm
;; compute driebes rescale legendre's


;; define factorial function

(define (factorial n)
   (define (iter product counter max-count)
      (if (> counter max-count)
         product
         (iter (* counter product)
                    (+ counter 1)
                    max-count )))
     (iter 1 1 n))

(define (pow x y)
    (exp (* y (log x))))

(define (I m n)
   (define (sgn l)  1)
   (define (term m n l)
           (/ (* (sgn l) (pow 0.5 l) (factorial (+ m n l)) )
            ( * (factorial (- (- n m) l)) 
                (factorial (+ (* 2 m) l 1)) 
                (factorial l) ) 
            )
           )
   (term m n 0))

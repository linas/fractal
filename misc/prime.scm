

 (define (prime n) (if (even? n) (= n 2) (let loop ((trial-divisor 3))
(cond ((< n (* trial-divisor trial-divisor)) (display "prime")) ((zero?
(remainder n trial-divisor))  (display trial-divisor)) (else (loop (+
trial-divisor 2)))))))




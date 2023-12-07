

; Markov numbers

(define mark (list 1 2 5 13 29 34 89 169 194 233 433 610 985 1325 1597
2897 4181 5741 6466 7561 9077 10946 14701 28657 33461 37666 43261 51641
62210 75025 96557 135137 195025))


(filter even? mark)
;; All odd markovs are 4k+1 for some k this returns k
;; k is not in OEIS
(map  (lambda (m) (/ (- m 1) 4)) (filter odd? mark))
; (0 1 3 7 22 42 58 108 246 331 399 724 1045 1435 1890 2269 3675 7164 8365 10815 12910 18756 24139 33784 48756)

;; All even markovs are 32j+2 some j
;; j is not in OEIS
(map  (lambda (m) (/ (- m 2) 32)) (filter even? mark))
; (0 1 6 19 202 342 1177 1944)

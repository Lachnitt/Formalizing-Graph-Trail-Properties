(declare-sort V 0)
(declare-sort W 0)

(declare-fun great (W W) Bool)
(declare-fun cost (V V) W)

(assert (forall ((x W)) (not (great x x))))
(assert (forall ((x W) (y W) (z W))
  (or (not (great x y)) (not (great y z)) (great x z))))
(assert (forall ((x W) (y W)) (or (great x y) (great y x) (= x y))))

(assert (forall ((x V) (y V)) (= (cost x y) (cost y x))))
(assert (forall ((w V) (x V) (y V) (z V))
  (or (and (= w x) (= y z)) (and (= y w) (= x z)) (not (= (cost x y) (cost w z))))))
(assert (forall ((x W)) (exists ((a V) (b V)) (= (cost a b) x))))


(assert (exists ((d0 V) (d1 V) (d2 V) (d3 V) (d4 V))
(and 
(distinct d0 d1 d2 d3 d4)
)))

(assert (forall ((c0 V) (c1 V) (c2 V) (c3 V) (c4 V) (c5 V))
  (not (and (great (cost c0 c1) (cost c1 c2))
            (distinct c0 c1)
            (distinct c1 c2)
	    (great (cost c1 c2) (cost c2 c3))
            (distinct c2 c3)
	    (great (cost c2 c3) (cost c3 c4))
            (distinct c3 c4)
	    (great (cost c3 c4) (cost c4 c5))
            (distinct c4 c5)
            ))))

(check-sat)

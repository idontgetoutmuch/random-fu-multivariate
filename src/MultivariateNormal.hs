module Data.Random.Distribution.MultivariateNormal
    ( Normal(..)
    ) where

normalMultivariate :: H.Vector Double -> H.Herm Double -> RVarT m (H.Vector Double)
normalMultivariate mu bigSigma = do
  z <- replicateM (H.size mu) (rvarT R.StdNormal)
  return $ mu + bigA H.#> (H.fromList z)
  where
    (vals, bigU) = H.eigSH bigSigma
    lSqrt = H.diag $ H.cmap sqrt vals
    bigA = bigU H.<> lSqrt

data family Normal k :: *

data instance Normal (H.Vector Double) = Normal (H.Vector Double) (H.Herm Double)

instance Distribution Normal (H.Vector Double) where
  rvar (Normal m s) = normalMultivariate m s

normalPdf :: (H.Numeric a, H.Field a, H.Indexable (H.Vector a) a, Num (H.Vector a)) =>
             H.Vector a -> H.Herm a -> H.Vector a -> a
normalPdf mu sigma x = exp $ normalLogPdf mu sigma x

normalLogPdf :: (H.Numeric a, H.Field a, H.Indexable (H.Vector a) a, Num (H.Vector a)) =>
                 H.Vector a -> H.Herm a -> H.Vector a -> a
normalLogPdf mu bigSigma x = - H.sumElements (H.cmap log (diagonals dec))
                              - 0.5 * (fromIntegral (H.size mu)) * log (2 * pi)
                              - 0.5 * s
  where
    dec = fromJust $ H.mbChol bigSigma
    t = fromJust $ H.linearSolve (H.tr dec) (H.asColumn $ x - mu)
    u = H.cmap (\x -> x * x) t
    s = H.sumElements u

diagonals :: (Storable a, H.Element t, H.Indexable (H.Vector t) a) =>
             H.Matrix t -> H.Vector a
diagonals m = H.fromList (map (\i -> m H.! i H.! i) [0..n-1])
  where
    n = max (H.rows m) (H.cols m)

instance PDF Normal (H.Vector Double) where
  pdf (Normal m s) = normalPdf m s
  logPdf (Normal m s) = normalLogPdf m s

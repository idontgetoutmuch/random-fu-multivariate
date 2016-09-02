-----------------------------------------------------------------------------
-- |
-- Module      :  Data.Random.Distribution.MultivariateNormal
-- Copyright   :  (c) 2016 FP Complete Corporation
-- License     :  MIT (see LICENSE)
-- Maintainer  :  dominic@steinitz.org
--
-- Sample from the multivariate normal distribution with a given
-- vector-valued \(\mu\) and covariance matrix \(\Sigma\). For example,
-- the chart below shows samples from the bivariate normal
-- distribution.
--
-- <<diagrams/src_Data_Random_Distribution_MultivariateNormal_diagM.svg#diagram=diagM&height=600&width=500>>
--
-- Example code to generate the chart:
--
-- > import qualified Graphics.Rendering.Chart as C
-- > import Graphics.Rendering.Chart.Backend.Diagrams
-- >
-- > import Data.Random.Distribution.MultivariateNormal
-- >
-- > import qualified Data.Random as R
-- > import Data.Random.Source.PureMT
-- > import Control.Monad.State
-- > import qualified Numeric.LinearAlgebra.HMatrix as LA
-- >
-- > nSamples :: Int
-- > nSamples = 10000
-- >
-- > sigma1, sigma2, rho :: Double
-- > sigma1 = 3.0
-- > sigma2 = 1.0
-- > rho = 0.5
-- >
-- > singleSample :: R.RVarT (State PureMT) (LA.Vector Double)
-- > singleSample = R.sample $ Normal (LA.fromList [0.0, 0.0])
-- >                (LA.sym $ (2 LA.>< 2) [ sigma1, rho * sigma1 * sigma2
-- >                                      , rho * sigma1 * sigma2, sigma2])
-- >
-- > multiSamples :: [LA.Vector Double]
-- > multiSamples = evalState (replicateM nSamples $ R.sample singleSample) (pureMT 3)
-- > pts = map (f . LA.toList) multiSamples
-- >   where
-- >     f [x, y] = (x, y)
-- >     f _      = error "Only pairs for this chart"
-- >
-- >
-- > chartPoint pointVals n = C.toRenderable layout
-- >   where
-- >
-- >     fitted = C.plot_points_values .~ pointVals
-- >               $ C.plot_points_style  . C.point_color .~ opaque red
-- >               $ C.plot_points_title .~ "Sample"
-- >               $ def
-- >
-- >     layout = C.layout_title .~ "Sampling Bivariate Normal (" ++ (show n) ++ " samples)"
-- >            $ C.layout_y_axis . C.laxis_generate .~ C.scaledAxis def (-3,3)
-- >            $ C.layout_x_axis . C.laxis_generate .~ C.scaledAxis def (-3,3)
-- >
-- >            $ C.layout_plots .~ [C.toPlot fitted]
-- >            $ def
-- >
-- > diagM = do
-- >   denv <- defaultEnv C.vectorAlignmentFns 600 500
-- >   return $ fst $ runBackend denv (C.render (chartPoint pts nSamples) (500, 500))
--
-----------------------------------------------------------------------------

{-# LANGUAGE TypeFamilies          #-}
{-# LANGUAGE FlexibleInstances     #-}
{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module Data.Random.Distribution.MultivariateNormal
    ( Normal(..)
    ) where

import           Data.Random.Distribution
import qualified Numeric.LinearAlgebra.HMatrix as H
import           Control.Monad
import qualified Data.Random as R
import           Foreign.Storable ( Storable )
import           Data.Maybe ( fromJust )

normalMultivariate :: H.Vector Double -> H.Herm Double -> R.RVarT m (H.Vector Double)
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
    u = H.cmap (\v -> v * v) t
    s = H.sumElements u

diagonals :: (Storable a, H.Element t, H.Indexable (H.Vector t) a) =>
             H.Matrix t -> H.Vector a
diagonals m = H.fromList (map (\i -> m H.! i H.! i) [0..n-1])
  where
    n = max (H.rows m) (H.cols m)

instance PDF Normal (H.Vector Double) where
  pdf (Normal m s) = normalPdf m s
  logPdf (Normal m s) = normalLogPdf m s

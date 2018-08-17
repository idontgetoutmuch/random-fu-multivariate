-----------------------------------------------------------------------------
-- |
-- Module      :  Data.Random.Distribution.Static.MultivariateNormal
-- Copyright   :  (c) 2016 FP Complete Corporation
-- License     :  MIT (see LICENSE)
-- Maintainer  :  dominic@steinitz.org
--
-- Sample from the multivariate normal distribution with a given
-- vector-valued \(\mu\) and covariance matrix \(\Sigma\). For
-- example, the chart below shows samples from the bivariate normal
-- distribution. The dimension of the mean \(n\) is statically checked
-- to be compatible with the dimension of the covariance matrix \(n \times n\).
--
-- <<diagrams/src_Data_Random_Distribution_Static_MultivariateNormal_diagMS.svg#diagram=diagMS&height=600&width=500>>
--
-- Example code to generate the chart:
--
-- > {-# LANGUAGE DataKinds #-}
-- >
-- > import qualified Graphics.Rendering.Chart as C
-- > import Graphics.Rendering.Chart.Backend.Diagrams
-- >
-- > import Data.Random.Distribution.Static.MultivariateNormal
-- >
-- > import qualified Data.Random as R
-- > import Data.Random.Source.PureMT
-- > import Control.Monad.State
-- > import Numeric.LinearAlgebra.Static
-- >
-- > nSamples :: Int
-- > nSamples = 10000
-- >
-- > sigma1, sigma2, rho :: Double
-- > sigma1 = 3.0
-- > sigma2 = 1.0
-- > rho = 0.5
-- >
-- > singleSample :: R.RVarT (State PureMT) (R 2)
-- > singleSample = R.sample $ Normal (vector [0.0, 0.0])
-- >                (sym $ matrix [ sigma1, rho * sigma1 * sigma2
-- >                              , rho * sigma1 * sigma2, sigma2])
-- >
-- > multiSamples :: [R 2]
-- > multiSamples = evalState (replicateM nSamples $ R.sample singleSample) (pureMT 3)
-- >
-- > pts = map f multiSamples
-- >   where
-- >     f z = (x, y)
-- >       where
-- >         (x, t) = headTail z
-- >         (y, _) = headTail t
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
-- > diagMS = do
-- >   denv <- defaultEnv C.vectorAlignmentFns 600 500
-- >   return $ fst $ runBackend denv (C.render (chartPoint pts nSamples) (500, 500))
--
-----------------------------------------------------------------------------

{-# OPTIONS_GHC -Wall                     #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing  #-}
{-# OPTIONS_GHC -fno-warn-type-defaults   #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind  #-}
{-# OPTIONS_GHC -fno-warn-missing-methods #-}
{-# OPTIONS_GHC -fno-warn-orphans         #-}

{-# LANGUAGE MultiParamTypeClasses        #-}
{-# LANGUAGE TypeFamilies                 #-}
{-# LANGUAGE ScopedTypeVariables          #-}
{-# LANGUAGE DataKinds                    #-}

module Data.Random.Distribution.Static.MultivariateNormal
    ( Normal(..)
    ) where

import           Data.Random hiding ( StdNormal, Normal )
import qualified Data.Random as R
import           Control.Monad.State ( replicateM )
import qualified Numeric.LinearAlgebra.HMatrix as H
import           Numeric.LinearAlgebra.Static as S
                 ( R, vector, extract, Sq, Sym, col,
                   tr, linSolve, uncol, chol, (<.>),
                   ℝ, (<>), diag, (#>), eigensystem
                 )
import          GHC.TypeLits ( KnownNat, natVal )
import          Data.Maybe ( fromJust )


normalMultivariate :: KnownNat n =>
                      R n -> Sym n -> RVarT m (R n)
normalMultivariate mu bigSigma = do
  z <- replicateM (fromIntegral $ natVal mu) (rvarT R.StdNormal)
  return $ mu + bigA #> (vector z)
  where
    (vals, bigU) = eigensystem bigSigma
    lSqrt = diag $ mapVector sqrt vals
    bigA = bigU S.<> lSqrt

mapVector :: KnownNat n => (ℝ -> ℝ) -> R n -> R n
mapVector f = vector . H.toList . H.cmap f . extract

sumVector :: KnownNat n => R n -> ℝ
sumVector x = x <.> 1

data family Normal k :: *

data instance Normal (R n) = Normal (R n) (Sym n)

instance KnownNat n => Distribution Normal (R n) where
  rvar (Normal m s) = normalMultivariate m s

normalLogPdf :: KnownNat n =>
                R n -> Sym n -> R n -> Double
normalLogPdf mu bigSigma x = - sumVector (mapVector log (diagonals dec))
                             - 0.5 * (fromIntegral $ natVal mu) * log (2 * pi)
                             - 0.5 * s
  where
    dec = chol bigSigma
    t = uncol $ fromJust $ linSolve (tr dec) (col $ x - mu)
    u = mapVector (\x -> x * x) t
    s = sumVector u

normalPdf :: KnownNat n =>
             R n -> Sym n -> R n -> Double
normalPdf mu sigma x = exp $ normalLogPdf mu sigma x

diagonals :: KnownNat n => Sq n -> R n
diagonals = vector . H.toList . H.takeDiag . extract

instance KnownNat n => PDF Normal (R n) where
  pdf (Normal m s) = normalPdf m s
  logPdf (Normal m s) = normalLogPdf m s

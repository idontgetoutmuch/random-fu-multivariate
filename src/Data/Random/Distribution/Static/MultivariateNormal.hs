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

module Data.Random.Distribution.MultivariateNormal.Static
    ( Normal(..)
    ) where

import           Data.Random hiding ( StdNormal, Normal )
import qualified Data.Random as R
import           Control.Monad.State ( replicateM )
import qualified Numeric.LinearAlgebra.HMatrix as H
import           Numeric.LinearAlgebra.Static
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
    bigA = bigU <> lSqrt

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

instance KnownNat n => PDF MultivariateNormal (R n) where
  pdf (MultivariateNormal m s) = normalPdf m s
  logPdf (MultivariateNormal m s) = normalLogPdf m s

--------------------------------------------------------------
--- An implementation of multivariate normal distributions ---
--------------------------------------------------------------
{-
Written by: Dominic Steinitz, Jacob West
Last modified: 2016-07-27

Summary: Multivariate normal distributions are necessary for Kalman
filters and smoothers.  However, strictly speaking, the functionality
provided here should exist elsewhere, perhaps in the package:
random-fu.
-}

---------------------------
--- File header pragmas ---
---------------------------
{-# LANGUAGE RecordWildCards #-}       -- Used by multiNormalRV, multiNormalConstant
                                       -- and multiNormalQuadraticForm
{-# LANGUAGE MultiParamTypeClasses #-} -- Necessary for Distribution instance
{-# LANGUAGE FlexibleInstances #-}     -- Necessary for Show instance
{-# LANGUAGE TypeFamilies #-}          -- Necessary for MultiNormal definition

------------------------
--- Module / Exports ---
------------------------
module Data.Random.Distribution.MultiNormal
       (
         MultiNormal(..)
       , inv
       )
       where
---------------
--- Imports ---
---------------
import Control.Monad (replicateM, when)
import Data.Maybe (fromMaybe, fromJust)
import Data.Random
import GHC.TypeLits
import Numeric.LinearAlgebra.Static

import qualified Numeric.LinearAlgebra as LA

------------------------
--- Helper Functions ---
------------------------
-- Matrix inverse: for some reason, this isn't built into the
-- static interface; warning: no error handling

-- WARNING: Needs better error handling
inv :: KnownNat n => Sq n -> Sq n
inv = fromMaybe (error "Failed attempting to invert non-invertible matrix.") .
      flip linSolve eye

----------------------------------------
--- Multivariate Normal Distrubtions ---
----------------------------------------
-- This probably belongs elsewhere, maybe Data.Random, but that would
-- create a dependence on Numric.LinearAlgebra which I believe is not
-- there now and may be undesirable.

data family MultiNormal k :: *
data instance KnownNat n => MultiNormal (R n) =
  MultiNormal { mu :: (R n), cov :: (Sym n) }

--- Show Instance ---
instance KnownNat n => Show (MultiNormal (R n)) where
  show MultiNormal {..} = "Normal " ++ show mu ++ " " ++ show cov

--- Distribution Instance ---
instance KnownNat n => Distribution MultiNormal (R n) where
  rvar = multiNormalRV

-- WARNING: Needs better error handling
multiNormalRV :: KnownNat n => MultiNormal (R n) -> RVarT m (R n)
multiNormalRV MultiNormal {..} = do
  let (vals, vecs) = eigensystem cov
  when (any (<0) (LA.toList $ unwrap vals))
    (error "Covariance matrix is not positive semi-definite.")

  let lSqrt = diag (fromJust . create $ LA.cmap sqrt (extract vals))
      bigA  = tr vecs <> lSqrt
  
  gnoise <- replicateM (size mu) (rvarT StdNormal)
  return $ mu + bigA #> (vector gnoise)

--- PDF Instance ---
instance KnownNat n => PDF MultiNormal (R n) where
  pdf    = multiNormalPDF
  logPdf = multiNormalLogPDF

multiNormalPDF :: KnownNat n => MultiNormal (R n) -> R n -> Double
multiNormalPDF mn pt =
  multiNormalConstant mn * exp (multiNormalQuadraticForm mn pt)

multiNormalLogPDF :: KnownNat n => MultiNormal (R n) -> R n -> Double
multiNormalLogPDF mn pt =
  multiNormalConstant mn + multiNormalQuadraticForm mn pt

multiNormalConstant :: KnownNat n => MultiNormal (R n) -> Double
multiNormalConstant MultiNormal {..} = recip . sqrt $ (2*pi)^n * detCov
  where
    n = size mu
    detCov = LA.det . extract . unSym $ cov

multiNormalQuadraticForm :: KnownNat n => MultiNormal (R n) -> R n -> Double
multiNormalQuadraticForm MultiNormal {..} pt = (diff LA.<.> invCov LA.#> diff) / (-2)
  where
    diff = extract (mu - pt)
    invCov = extract . inv . unSym $ cov

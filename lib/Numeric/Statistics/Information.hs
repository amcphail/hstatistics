{-# LANGUAGE FlexibleContexts #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Numeric.Statistics.Information
-- Copyright   :  (c) A. V. H. McPhail 2010, 2014
-- License     :  BSD3
--
-- Maintainer  :  haskell.vivian.mcphail <at> gmail <dot> com
-- Stability   :  provisional
-- Portability :  portable
--
-- Shannon entropy
--
-----------------------------------------------------------------------------

module Numeric.Statistics.Information (
                                   entropy
                                   , mutual_information
                                  ) where


-----------------------------------------------------------------------------

import Numeric.Statistics.PDF

import qualified Numeric.LinearAlgebra as LA
--import Numeric.LinearAlgebra.Data hiding(Vector)

--import qualified Data.Vector as DV
import qualified Data.Vector.Storable as V

import Prelude hiding(map,zip)

-----------------------------------------------------------------------------

zeroToOne x
    | x == 0.0  = 1.0
    | otherwise = x

logE = V.map (log . zeroToOne)


-----------------------------------------------------------------------------

-- | the entropy \sum p_i l\ln{p_i} of a sequence
entropy :: PDF a Double 
        => a                       -- ^ the underlying distribution
        -> LA.Vector Double           -- ^ the sequence
        -> Double                  -- ^ the entropy
entropy p x = let ps = probability p x
              in negate $ (LA.dot ps (logE ps))

-- | the mutual information \sum_x \sum_y p(x,y) \ln{\frac{p(x,y)}{p(x)p(y)}}
mutual_information :: (PDF a Double, PDF b (Double,Double)) 
                   => b                                          -- ^ the underlying distribution
                   -> a                                          -- ^ the first dimension distribution
                   -> a                                          -- ^ the second dimension distribution
                   -> (LA.Vector Double, LA.Vector Double)             -- ^ the sequence
                   -> Double         -- ^ the mutual information
mutual_information p px py (x,y) = let ps = probability p $ V.zipWith (,) x y
                                       xs = probability px x
                                       ys = probability py y
                                   in (LA.dot ps (logE ps - logE (xs*ys)))

-----------------------------------------------------------------------------

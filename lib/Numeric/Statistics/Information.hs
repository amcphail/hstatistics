{-# OPTIONS_GHC -fglasgow-exts #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Numeric.Statistics.Information
-- Copyright   :  (c) A. V. H. McPhail 2010
-- License     :  GPL-style
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

import Data.Packed.Vector

import Numeric.Statistics.PDF

import Numeric.LinearAlgebra.Algorithms()
import Numeric.LinearAlgebra.Interface()
import Numeric.LinearAlgebra.Linear

-----------------------------------------------------------------------------

zeroToOne x
    | x == 0.0  = 1.0
    | otherwise = x

logE = mapVector (log . zeroToOne)


-----------------------------------------------------------------------------

-- | the entropy \sum p_i l\ln{p_i} of a sequence
entropy :: PDF a Double 
        => a                       -- ^ the underlying distribution
        -> Vector Double           -- ^ the sequence
        -> Double                  -- ^ the entropy
entropy p x = let ps = probability p x
              in negate $ dot ps (logE ps)

-- | the mutual information \sum_x \sum_y p(x,y) \ln{\frac{p(x,y)}{p(x)p(y)}}
mutual_information :: (PDF a Double, PDF b (Double,Double)) 
                   => b                                          -- ^ the underlying distribution
                   -> a                                          -- ^ the first dimension distribution
                   -> a                                          -- ^ the second dimension distribution
                   -> (Vector Double, Vector Double)             -- ^ the sequence
                   -> Double         -- ^ the mutual information
mutual_information p px py (x,y) = let ps = probability p $ zipVector x y
                                       xs = probability px x
                                       ys = probability py y
                                   in negate $ dot ps (logE ps - logE (xs*ys)) 

-----------------------------------------------------------------------------

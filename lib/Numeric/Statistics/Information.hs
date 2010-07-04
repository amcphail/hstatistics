{-# OPTIONS_GHC -fglasgow-exts #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Numeric.Statistics.Information
-- Copyright   :  (c) Alexander Vivian Hugh McPhail 2010
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

import qualified Numeric.GSL.Histogram as H
import qualified Numeric.GSL.Histogram2D as H2

import Numeric.LinearAlgebra.Algorithms
import Numeric.LinearAlgebra.Interface()

-----------------------------------------------------------------------------

zeroToOne x
    | x == 0.0  = 1.0
    | otherwise = x

logE = mapVector (log . zeroToOne)


-----------------------------------------------------------------------------

-- | the entropy \sum p_i l\ln{p_i} of a sequence
entropy :: H.Histogram             -- ^ the underlying distribution
        -> Vector Double           -- ^ the sequence (expected to fall within bounds of Histogram)
        -> Double                  -- ^ the entropy
entropy p x = let ps = H.prob p x
              in negate $ dot ps (logE ps)

-- | the mutual information \sum_x \sum_y \ln{\frac{p(x,y)}{p(x)p(y)}}
mutual_information :: H2.Histogram2D -- ^ the underlying distribution
                   -> H.Histogram    -- ^ the first dimension distribution
                   -> H.Histogram    -- ^ the second dimension distribution
                   -> (Vector Double, Vector Double) -- ^ the sequence (expected to fall within bounds of Histogram)
                   -> Double         -- ^ the mutual information
mutual_information p px py z@(x,y) = let ps = H2.prob p z
                                         xs = H.prob px x
                                         ys = H.prob py y
                               in negate $ dot ps (logE ps - logE (xs*ys)) 

-----------------------------------------------------------------------------

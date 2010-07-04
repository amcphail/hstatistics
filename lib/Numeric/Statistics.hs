{-# OPTIONS_GHC -fglasgow-exts #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Numeric.Statistics
-- Copyright   :  (c) Alexander Vivian Hugh McPhail 2010
-- License     :  GPL-style
--
-- Maintainer  :  haskell.vivian.mcphail <at> gmail <dot> com
-- Stability   :  provisional
-- Portability :  portable
--
-- Useful statistical functions
--
-----------------------------------------------------------------------------

module Numeric.Statistics (
                          covarianceMatrix
                          ) where


-----------------------------------------------------------------------------

import Data.Packed.Vector
import Data.Packed.Matrix

import qualified Data.Array.IArray as I 

import Numeric.GSL.Statistics

-----------------------------------------------------------------------------

-- | the covariance matrix
covarianceMatrix :: I.Array Int (Vector Double) -- ^ the dimensions of data (each vector being one dimension)
                 -> Matrix Double               -- ^ the symmetric covariance matrix
covarianceMatrix d = let (s,f) = I.bounds d
                      in fromArray2D $ I.array ((s,s),(f,f)) $ concat $ map (\(x,y) -> let c = covariance (d I.! x) (d I.! y) in if x == y then [((x,y),c)] else [((x,y),c),((y,x),c)]) $ filter (\(x,y) -> x <= y) $ I.range ((s,s),(f,f))

-----------------------------------------------------------------------------

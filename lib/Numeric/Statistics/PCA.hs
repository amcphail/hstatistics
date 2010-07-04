{-# OPTIONS_GHC -fglasgow-exts #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Numeric.Statistics.PCA
-- Copyright   :  (c) Alexander Vivian Hugh McPhail 2010
-- License     :  GPL-style
--
-- Maintainer  :  haskell.vivian.mcphail <at> gmail <dot> com
-- Stability   :  provisional
-- Portability :  portable
--
-- Principle Components Analysis
--
-----------------------------------------------------------------------------

module Numeric.Statistics.PCA (
                               pca, pcaTransform
                          ) where


-----------------------------------------------------------------------------

import qualified Data.Array.IArray as I 

import Data.Packed.Vector
import Data.Packed.Matrix

import Numeric.LinearAlgebra.Interface
import Numeric.LinearAlgebra.Algorithms

import Numeric.GSL.Statistics

import Numeric.Statistics

-----------------------------------------------------------------------------

qcomp e1@(ev1,_) e2@(ev2,_) 
    | ev1 < ev2             = True
    | otherwise             = False

qsort (x:xs) = [ y | y <- xs, qcomp y x ] ++ (x : [ y | y <- xs, qcomp x y ])

-- | find the n principal components of multidimensional data
pca :: I.Array Int (Vector Double)    -- the data
    -> Int                            -- how many dimensions to keep
    -> Matrix Double
pca d n = let d' = fmap (\x -> x - (scalar $ mean x)) d -- remove the mean from each dimension
              cv = covarianceMatrix d'
              (val',vec') = eigSH cv           -- the covariance matrix is real symmetric
              val = toList val'
              vec = toColumns vec'
              v' = zip val vec
              v = take n $ qsort v'
          in fromColumns $ snd $ unzip v

-- | perform a PCA transform of the original data (remove mean)
-- |     Final = M^T Data^T
pcaTransform :: I.Array Int (Vector Double)    -- ^ the data
             -> Matrix Double                  -- ^ the principal components
             -> I.Array Int (Vector Double)    -- ^ the transformed data
pcaTransform d m = let d' = fmap (\x -> x - (scalar $ mean x)) d -- remove the mean from each dimension
                       a = fromRows $ I.elems d'
                   in I.listArray (1,cols m) $ toColumns $ (trans m) <> a

-----------------------------------------------------------------------------

-----------------------------------------------------------------------------
-- |
-- Module      :  Numeric.Statistics.Surrogate
-- Copyright   :  (c) Alexander Vivian Hugh McPhail 2010
-- License     :  BSD3
--
-- Maintainer  :  haskell.vivian.mcphail <at> gmail <dot> com
-- Stability   :  provisional
-- Portability :  portable
--
-- Methods for tests using surrogate data
--
-----------------------------------------------------------------------------

module Numeric.Statistics.Surrogate (
                                     surrogate
                                    ) where


-----------------------------------------------------------------------------

import Data.Packed.Vector
--import Data.Packed.Matrix

import qualified Data.Array.IArray as I 

import Numeric.GSL.Permutation

import System.Random

-----------------------------------------------------------------------------

-- | perform an analysis using surrogate data
surrogate :: Int                          -- ^ random seed
          -> Int                          -- ^ number of repetitions
          -> (I.Array Int (Vector Double) -> a) -- ^ the evaluation function
          -> I.Array Int (Vector Double)  -- ^ the data
          -> I.Array Int a                -- ^ the results, with the evaluated real data in position 1 and the rest of the array containing the evaluated surrogate data
                                                                                                              
surrogate r n f d = I.listArray (1,n+1) $ (f d) : (surrogate' (mkStdGen r) n f d)

surrogate' :: StdGen -> Int-> (I.Array Int (Vector Double) -> a) -> I.Array Int (Vector Double) -> [a]
surrogate' _ 0 _ _ = []
surrogate' g n f d = let (g',g'') = split g
                         d' = permute_data g' d
                     in (f d) : (surrogate' g'' (n-1) f d')

randomList :: StdGen -> Int -> [Int]
randomList _ 0 = []
randomList g n = let (r,g') = random g
                 in r : (randomList g' (n-1))

permute_data :: StdGen -> I.Array Int (Vector Double) -> I.Array Int (Vector Double)
permute_data g d = let s = I.rangeSize $ I.bounds d
                       rs = randomList g s
                       ds = zip rs $ I.elems d
                   in I.listArray (I.bounds d) $ map (\(r,v) -> permute (random_permute r (dim v)) v) ds
 
-----------------------------------------------------------------------------

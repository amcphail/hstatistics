{-# LANGUAGE FlexibleContexts #-}
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
                           Sample,Samples
                          , covarianceMatrix
                          , meanList, meanArray, meanMatrix
                          , varianceList, varianceArray, varianceMatrix
                          ) where


-----------------------------------------------------------------------------

--import Numeric.Vector
--import Numeric.Matrix
--import Numeric.Container
import Numeric.LinearAlgebra

import qualified Data.Array.IArray as I 

import Numeric.GSL.Statistics

-----------------------------------------------------------------------------

type Sample a = Vector a
type Samples a = I.Array Int (Vector a)

-----------------------------------------------------------------------------

-- | the covariance matrix
covarianceMatrix :: Samples Double              -- ^ the dimensions of data (each vector being one dimension)
                 -> Matrix Double               -- ^ the symmetric covariance matrix
covarianceMatrix d = let (s,f) = I.bounds d
                      in fromArray2D $ I.array ((s,s),(f,f)) $ concat $ map (\(x,y) -> let c = covariance (d I.! x) (d I.! y) in if x == y then [((x,y),c)] else [((x,y),c),((y,x),c)]) $ filter (\(x,y) -> x <= y) $ I.range ((s,s),(f,f))

-----------------------------------------------------------------------------

-- | the mean of a list of vectors
meanList :: (Container Vector a, Num (Vector a)) => [Sample a] -> Sample a
meanList []     = error "meanVectors: empty list"
meanList [s]    = s
meanList (s:ss) = let ln = fromIntegral $ length ss + 1
                  in scale (recip ln) $ foldl (+) s ss

-- | the mean of an array of vectors
meanArray :: (Container Vector a, Num (Vector a)) => Samples a -> Sample a
meanArray a = meanList $ I.elems a

-- | the mean of a matrix with data series in rows
meanMatrix :: (Container Vector a, Num (Vector a), Element a) => Matrix a -> Sample a
meanMatrix a = meanList $ toRows a

-----------------------------------------------------------------------------

-- | the variance of a list of vectors
varianceList :: (Container Vector a, Floating (Vector a)) => [Sample a] -> Sample a
varianceList []  = error "varianceList: empty list"
varianceList [s] = constant 0 (dim s)
varianceList l   = let mxs = meanList (map (** 2) l)
                       msx = (meanList l) ** 2
                   in mxs - msx

-- | the variance of an array of vectors
varianceArray :: (Container Vector a, Floating (Vector a)) => Samples a -> Sample a
varianceArray a = varianceList $ I.elems a

-- | the variance of a matrix with data series in rows
varianceMatrix :: (Container Vector a, Floating (Vector a), Element a) => Matrix a -> Sample a
varianceMatrix a = varianceList $ toRows a

-----------------------------------------------------------------------------

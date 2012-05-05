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
                          , covarianceMatrix, correlationCoefficientMatrix
                          , meanList, meanArray, meanMatrix
                          , varianceList, varianceArray, varianceMatrix
                          --
                          , centre, cloglog, cut, ranks, kendall, logit
                          , mahanalobis
                          ) where


-----------------------------------------------------------------------------

--import Numeric.Vector
--import Numeric.Matrix
--import Numeric.Container
import Numeric.LinearAlgebra hiding (rank)

import qualified Data.Array.IArray as I 

import qualified Data.Vector.Generic as GV

import Numeric.GSL.Statistics
import Numeric.GSL.Sort(sort)

-----------------------------------------------------------------------------

type Sample a = Vector a
type Samples a = I.Array Int (Vector a)

-----------------------------------------------------------------------------

-- | the covariance matrix
covarianceMatrix :: Samples Double              -- ^ the dimensions of data (each vector being one dimension)
                 -> Matrix Double               -- ^ the symmetric covariance matrix
covarianceMatrix d = let (s,f) = I.bounds d
                      in fromArray2D $ I.array ((s,s),(f,f)) $ concat $ map (\(x,y) -> let c = covariance (d I.! x) (d I.! y) in if x == y then [((x,y),c)] else [((x,y),c),((y,x),c)]) $ filter (\(x,y) -> x <= y) $ I.range ((s,s),(f,f))

-- | the correlation coefficient: (cov x y) / (std x) (std y)
correlationCoefficientMatrix :: Samples Double -> Matrix Double
correlationCoefficientMatrix d = let (s,f) = I.bounds d
                           in fromArray2D $ I.array ((s,s),(f,f)) $ concat $ map (\(x,y) -> let { x' = d I.! x ; y' = d I.! y ; c = (covariance x' y') / ((stddev x') * (stddev y')) } in if x == y then [((x,y),c)] else [((x,y),c),((y,x),c)]) $ filter (\(x,y) -> x <= y) $ I.range ((s,s),(f,f))

-----------------------------------------------------------------------------
-----------------------------------------------------------------------------
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
-----------------------------------------------------------------------------
-----------------------------------------------------------------------------

-- | centre the data to 0: (x - (mean x))
centre :: Vector Double -> Vector Double
centre v = v - (realToFrac (mean v))

-----------------------------------------------------------------------------

-- | complementary log-log function
cloglog :: Vector Double -> Vector Double
cloglog v = - log (- (log v))

-----------------------------------------------------------------------------

-- | cut numerical data into intervals, data must fall inside the bounds
cut :: Vector Double 
    -> Vector Double -- ^ intervals
    -> Vector Int    -- ^ data indexed by bin
cut v c  = let c' = sort c
           in mapVector (\x -> cut_helper 0 x c') v 
    where
      cut_helper i x c 
          | i >= dim c                       = error "Numeric.Statistics: cut: data point not within interval"
          | x >= (c @> i) && x <= (c @> (i+1)) = i
          | otherwise                       = cut_helper (i + 1) x c

-----------------------------------------------------------------------------

-- | return the rank of each element of the vector
--     multiple identical entries result in the average rank of those entries
ranks :: Vector Double -> Vector Double
ranks v = let v' = sort v
          in mapVector (\x -> 1 + rank_helper x v') v
              where rank_helper x v' = let is = GV.elemIndices x v'
                                       in (realToFrac (GV.foldl (+) 0 is)) / (fromIntegral $ dim is)

-----------------------------------------------------------------------------

-- | kendall's rank correlation τ
kendall :: Vector Double -> Vector Double -> Matrix Double
kendall x y = let ln = dim x
                  rx = ranks x
                  ry = ranks y
                  r = fromColumns [rx,ry]
                  m = signum $ (kronecker r (asColumn $ constant 1.0 ln)) - (kronecker (asRow $ constant 1.0 ln) r)
                  c = rows m - 1
              in correlationCoefficientMatrix $ I.listArray (0,c) (toColumns m)

-----------------------------------------------------------------------------

-- | (logit p) = log(p/(1-p))
logit :: Vector Double -> Vector Double
logit v =  mapVector (\x -> - (log ((1 / x) - 1))) v

-----------------------------------------------------------------------------

-- | the Mahalanobis D-square distance between samples
--     columns are components and rows are observations
mahalanobis :: Samples Double        -- ^ the data set
            -> Maybe (Sample Double) -- ^ (Just sample) to be measured or use mean when Nothing
            -> Double                -- ^ D^2 
mahalanobis x u = let (_,xr) = I.bounds x
                      xl     = I.elems x
                      s'     = pinv $ covarianceMatrix x
                      xu     = case u of
                                 Nothing -> fromList $ map mean xl
                                 Just m  -> m
                      xm     = fromRows $ map ((-) xu) $ toRows $ fromColumns xl
                      --um     = asColumn xu
                      --w      = ((trans xm) <> xm + (trans um) <> um)/(fromIntegral $ xr - 1)
                      --w'     = inv w
                  in ((xm <> s' <> (trans xm)) @@> (0,0)) 

m_test_samples :: Samples Double
m_test_samples = I.listArray (1,3) $ map fromList [[1,2,3,4],[2,3,4,5],[1,3,2,4]]

-----------------------------------------------------------------------------


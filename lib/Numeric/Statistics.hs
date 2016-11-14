{-# LANGUAGE FlexibleContexts #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Numeric.Statistics
-- Copyright   :  (c) A. V. H. McPhail 2010, 2012, 2014
-- License     :  BSD3
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
                          , centre, cloglog, corcoeff, cut
                          , ranks, kendall, logit
                          , mahalanobis
                          , mode, moment
                          , ols, percentile, range
                          , run_count
                          , spearman, studentize
                          ) where


-----------------------------------------------------------------------------

import Numeric.LinearAlgebra hiding(range)
--import Numeric.LinearAlgebra.Data hiding(range)
--import Numeric.LinearAlgebra.Devel

import qualified Data.Array.IArray as I 
import qualified Data.List as DL
import qualified Data.Vector.Generic as GV
--import qualified Data.Vector.Storable as SV

import Foreign.Storable

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
meanList :: (Container Vector a, Num (Vector a), Fractional a) => [Sample a] -> Sample a
meanList []     = error "meanVectors: empty list"
meanList [s]    = s
meanList (s:ss) = let ln = fromIntegral $ length ss + 1
                  in scale (recip ln) $ foldl (+) s ss

-- | the mean of an array of vectors
meanArray :: (Container Vector a, Num (Vector a), Fractional a) => Samples a -> Sample a
meanArray a = meanList $ I.elems a

-- | the mean of a matrix with data series in rows
meanMatrix :: (Container Vector a, Num (Vector a), Fractional a) => Matrix a -> Sample a
meanMatrix a = meanList $ toRows a

-----------------------------------------------------------------------------

-- | the variance of a list of vectors
varianceList :: (Container Vector a, Floating (Vector a), Num a, Fractional a) => [Sample a] -> Sample a
varianceList []  = error "varianceList: empty list"
varianceList [s] = konst 0 (size s)
varianceList l   = let mxs = meanList (map (** 2) l)
                       msx = (meanList l) ** 2
                   in mxs - msx

-- | the variance of an array of vectors
varianceArray :: (Container Vector a, Floating (Vector a), Fractional a) => Samples a -> Sample a
varianceArray a = varianceList $ I.elems a

-- | the variance of a matrix with data series in rows
varianceMatrix :: (Container Vector a, Floating (Vector a), Fractional a) => Matrix a -> Sample a
varianceMatrix a = varianceList $ toRows a

-----------------------------------------------------------------------------
-----------------------------------------------------------------------------
-----------------------------------------------------------------------------

-- | centre the data to 0: (x - (mean x))
centre :: Vector Double -> Vector Double
centre v = v - (realToFrac (mean v))

-----------------------------------------------------------------------------

-- | complementary log-log function
--cloglog :: Vector Double -> Vector Double
cloglog :: Floating a => a -> a
cloglog v = - log (- (log v))

-----------------------------------------------------------------------------

-- | corcoeff = covariance x / (std dev x * std dev y)
corcoeff :: Vector Double -> Vector Double -> Double
corcoeff x y = (covariance x y)/((stddev x)*(stddev y))

-----------------------------------------------------------------------------

-- | cut numerical data into intervals, data must fall inside the bounds
cut :: Vector Double 
    -> Vector Double -- ^ intervals
    -> Vector Int    -- ^ data indexed by bin
cut v c  = let c' = sort c
           in GV.map (\x -> cut_helper 0 x c') v 
    where
      cut_helper j x d 
          | j >= size d                       = error "Numeric.Statistics: cut: data point not within interval"
          | x >= (d `atIndex` j) && x <= (d `atIndex` (j+1)) = j
          | otherwise                       = cut_helper (j + 1) x d

-----------------------------------------------------------------------------

-- | return the rank of each element of the vector
--     multiple identical entries result in the average rank of those entries
--ranks :: Vector Double -> Vector Double
ranks :: (Fractional b, Storable b) => Vector Double -> Vector b
ranks v = let v' = sort v
          in GV.map (\x -> 1 + rank_helper x v') v
              where rank_helper x v' = let is = GV.elemIndices x v'
                                       in (realToFrac (GV.foldl (+) 0 is)) / (fromIntegral $ GV.length is)

-----------------------------------------------------------------------------

-- | kendall's rank correlation τ
kendall :: Vector Double -> Vector Double -> Matrix Double
kendall x y = let ln = size x
                  rx = ranks x
                  ry = ranks y
                  r = fromColumns [rx,ry]
                  m = signum $ (kronecker r (asColumn $ konst 1.0 ln)) - (kronecker (asRow $ konst 1.0 ln) r)
                  c = rows m - 1
              in correlationCoefficientMatrix $ I.listArray (0,c) (toColumns m)

-----------------------------------------------------------------------------

-- | (logit p) = log(p/(1-p))
--logit :: Vector Double -> Vector Double
logit :: (Floating b, Storable b)
        => Vector b -> Vector b
logit v =  GV.map (\x -> - (log ((1 / x) - 1))) v

-----------------------------------------------------------------------------

-- | the Mahalanobis D-square distance between samples
--     columns are components and rows are observations (uses pseudoinverse)
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
                      --w      = ((tr' xm) <> xm + (tr' um) <> um)/(fromIntegral $ xr - 1)
                      --w'     = inv w
                  in ((xm <> s' <> (tr' xm)) `atIndex` (0,0)) 

-----------------------------------------------------------------------------

-- | a list of element frequencies
mode :: Vector Double -> [(Double,Integer)]
mode v = let w = sort v
         in DL.sortBy (\(_,n) (_,n') -> compare n' n) $ GV.foldr freqs [] w
            where freqs x []          = [(x,1)]
                  freqs x ((f,n):fns)
                      | f == x         = ((f,n+1):fns) 
                      | otherwise     = ((x,1):(f,n):fns)

-----------------------------------------------------------------------------

-- | the p'th moment of a vector
moment :: Integral a 
       => a             -- ^ moment
       -> Bool          -- ^ calculate central moment
       -> Bool          -- ^ calculate absolute moment
       -> Vector Double -- ^ data
       -> Double
moment p c a v 
    | p <= 0     = error "Numeric.Statistics.moment: negative moment requested"
--    | p == 1     = mean v    
--    | p == 2     = variance v -- gives sample variance
    | otherwise = let u = if c then centre v else v
                      w = if a then abs u else u
                      x = GV.map (** (fromIntegral p)) w
                  in mean x

-----------------------------------------------------------------------------

-- | ordinary least squares estimation for the multivariate model
--   Y = X B + e        rows are observations, columns are elements
--   mean e = 0, cov e = kronecker s I
ols :: (Num (Vector t), Field t) 
      => Matrix t         -- ^ X
    -> Matrix t           -- ^ Y
    -> (Matrix t, Matrix t, Matrix t) -- ^ (OLS estimator for B, OLS estimator for s, OLS residuals)
ols x y 
    | rows x /= rows y = error "Numeric.Statistics: ols: incorrect matrix dimensions"
    | otherwise       = let (xr,xc) = (rows x,cols x)
                            (yr,yc) = (rows y,cols y)
                            z = (tr' x) <> x
                            r = rank z
                            beta = if r == xc 
                                      then (inv z) <> (tr' x) <> y
                                      else (pinv x) <> y
                            rr = y - x <> beta
                            sigma = ((tr' rr) <> rr) / (fromIntegral $ xr - r)
                        in (beta,rr,sigma)

-----------------------------------------------------------------------------

-- | compute quantiles in percent
percentile :: Double        -- ^ percentile (0 - 100)
           -> Vector Double -- ^ data
           -> Double        -- ^ result
percentile p d = quantile (0.01*p) d

-----------------------------------------------------------------------------

-- | the difference between the maximum and minimum of the input
range :: (Container c e, Num e) => c e -> e
range v = maxElement v - minElement v

-----------------------------------------------------------------------------

-- | count the number of runs greater than or equal to @n@ in the data
run_count :: (Num a, Num t, Ord b, Ord a, Container Vector b) 
            => a             -- ^ longest run to count
          -> Vector b        -- ^ data
          -> [(a, t)]        -- ^ [(run length,count)]
run_count n v = let w = subVector 1 (size v - 1) v
                    x = GV.foldr run_count' [(1,v `atIndex` 0)] w
                    y = map fst x
                    z = takeWhile (<= n) $  DL.sort y
                in foldr count [] z
    where run_count' m ((c,g):cs)
              | m < g             = ((c+1,m):cs)
              | otherwise         = ((1,m):(c,g):cs)
          count x []           = [(x,1)]
          count x ((yv,yc):ys)   
              | x == yv         = ((yv,yc+1):ys)
              | otherwise      = ((x,1):(yv,yc):ys)

-----------------------------------------------------------------------------

-- | Spearman's rank correlation coefficient
spearman :: Vector Double -> Vector Double -> Double
spearman x y = corcoeff (ranks x) (ranks y)

-----------------------------------------------------------------------------

-- | centre and normalise a vector
studentize :: Vector Double -> Vector Double
studentize x = (centre x)/(fromList $ [stddev x])

-----------------------------------------------------------------------------

--table

-----------------------------------------------------------------------------




-----------------------------------------------------------------------------


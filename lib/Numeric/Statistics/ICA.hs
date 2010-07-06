{-# OPTIONS_GHC -fglasgow-exts #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Numeric.Statistics.ICA
-- Copyright   :  (c) Alexander Vivian Hugh McPhail 2010
-- License     :  GPL-style
--
-- Maintainer  :  haskell.vivian.mcphail <at> gmail <dot> com
-- Stability   :  provisional
-- Portability :  portable
--
-- Independent Components Analysis
--
--  implements the FastICA algorithm found in:
--
--   http://www.google.com/url?sa=t&source=web&cd=2&ved=0CBgQFjAB&url=http%3A%2F%2Fciteseerx.ist.psu.edu%2Fviewdoc%2Fdownload%3Fdoi%3D10.1.1.79.7003%26rep%3Drep1%26type%3Dpdf&ei=RQozTJb6L4_fcbCV6cMD&usg=AFQjCNGClLIB9MAvbrEj45SyUx9cYubLyA&sig2=hg5Wnfy3dLPkoIc1hqSfjg
--
--   Aapo Hyvärinen and Erkki Oja
--   Independent Component Analysis: Algorithms and Applications
--   Neural Networks, 13(4-5):411-430, 2000
--
-----------------------------------------------------------------------------

module Numeric.Statistics.ICA (
                               sigmoid, sigmoid',
                               demean, whiten,
                               ica, icaDefaults
                          ) where


-----------------------------------------------------------------------------

import qualified Data.Array.IArray as I 

import Data.Packed.Vector
import Data.Packed.Matrix
--import Data.Packed.Random

import Numeric.LinearAlgebra.Interface
import Numeric.LinearAlgebra.Algorithms

import Numeric.GSL.Statistics

import Numeric.Statistics

import System.Random

-----------------------------------------------------------------------------

-- | sigmoid transfer function
sigmoid :: Double -> Double
sigmoid u = u * exp((-u**2)/2)

-- | derivative of sigmoid transfer function
sigmoid' :: Double -> Double
sigmoid' u = -u**2 * exp((-u**2)/2)

-----------------------------------------------------------------------------

-- preprocessing:
--   demean
--   whiten
--        eigenvalue decomposition of covariance matrix E{xx^T} = EDE^T
--                   E orthogonal matrix of eigenvectors
--                   D diagonal matrix of eigenvalues, D = diag(d_1,...,d_n)
--                   x_white = ED^{-1/2}E^Tx
--                   D^{-1/2} = diag{d_1^{-1/2},...}
--

-----------------------------------------------------------------------------

-- | remove the mean from data
demean :: I.Array Int (Vector Double)                  -- ^ the data
       -> (I.Array Int (Vector Double),Vector Double)  -- ^ (demeaned data,mean)
demean d = let u = I.elems $ fmap mean d
               d' = I.listArray (I.bounds d) (zipWith (-) (I.elems d) (map scalar u))
               u' = fromList u
               in (d',u')

-- | whiten data
whiten :: I.Array Int (Vector Double)                 -- ^ the data
       -> Double                                      -- ^ eigenvalue threshold
       -> (I.Array Int (Vector Double),Matrix Double) -- ^ (whitened data,transform)
whiten d q = let cv = covarianceMatrix d
                 (val',vec') = eigSH cv           -- the covariance matrix is real symmetric
                 val = toList val'
                 vec = toColumns vec'
                 v' = zip val vec
                 v = filter ((> q) . fst) v'        -- keep only eigens > than parameter
                 (dd',e') = unzip v
                 dd = diag $ (** (-0.5)) $ fromList dd'  -- square root of eigenvalues diagonalised
                 e = fromColumns e'
                 x = fromRows $ I.elems d
                 t = e <> dd <> trans e          -- the actual mathematics
                 x' = t <> x                     -- the actual mathematics
                 d' = I.listArray (I.bounds d) (toRows x') 
             in (d',t)
                 
-----------------------------------------------------------------------------

-- assuming that a weight vector is a row

-- algorithm:
-- 1  initial random weight vectors w_i
-- 2  w_i^+ = E{xg(w^Tx)} - E{g'(w^Tx)}w     (newton phase)
-- 3  W = W = (WW^T)^{-1/2)W                 (decorrelation) W = ( ..., w_i, ...)^T
--                                           WW^T = FDF^T (eigenvalue decomposition)
-- 4  w_i = w^+/norm(w^+)                    (normalisation) (almost any norm but not Frobenius)
-- 5  if not converged (dot w w^+ ~ 1 implies convergence) go to step 2
--
-- in matrix form, 2 becomes:
--    W^+ = W + (diag a_i)[(diag b_i) + E{g(y)y^T}]W
--
--      where
--            y = Wx
--            b_i = -E{y_ig(y_i)}
--            a_i = -1/(b_i-E{g'(y_i)})
--
--    g(u) = tanh(au) 0<=a<=2, often a = 1
--    g(u) = u exp(-u^2/2)

-----------------------------------------------------------------------------

unconcat 0     _ _  = []
unconcat (r+1) c xs = [take c xs] ++ unconcat r c (drop c xs)

random_vector :: Int -> (Int,Int) -> Matrix Double
random_vector s (r,c) = fromLists $ unconcat r c $ randomRs (-1,1) (mkStdGen s)

-- g g' w x -> w'
update :: (Double -> Double) -> (Double -> Double) -> Matrix Double -> Matrix Double -> Matrix Double
update g g' w x = let y = w <> x
                      ys = toRows y
                      bis = map (\y' -> - mean (y' * (mapVector g y'))) ys
                      ais = zipWith (\b y' -> -1 / (b - mean (mapVector g y'))) bis ys
                      r = rows y
                      ix = ((1,1),(r,r))
                      cov = fromArray2D $ I.listArray ix $ map (\(m,n) -> covariance (mapVector g' (ys!!(m-1))) (ys!!(n-1))) $ I.range ix
                  in w + (diag $ fromList ais) <> ((diag $ fromList bis) + cov) <> w  

decorrelate :: Matrix Double -> Matrix Double
decorrelate w = let w' = w / (scalar $ sqrt $ pnorm PNorm2 (w <> trans w))
                in decorrelate' w w'
    where decorrelate' w w' 
              | converged 0.000001 w w' = w'
              | otherwise               = decorrelate' w' ((scale 1.5 w') - (scale 0.5 (w <> trans w <> w)))
{- don't know how to do svd of non-square matrices
decorrelate m = let (u,d,v) = svd m
                in u <> (diag (d ** (-0.5))) <> trans v <> m
-}

normalise :: NormType -> Matrix Double -> Matrix Double
normalise t m = fromRows $ map (\v -> v / (scalar $ pnorm t v)) (toRows m)

converged :: Double -> Matrix Double -> Matrix Double -> Bool
converged t m m' = let d' = map ((-) 1) $ zipWith dot (toRows m) (toRows m')
                   in maximum d' <= t

-----------------------------------------------------------------------------

ica' :: (Double -> Double)          -- ^ transfer function (tanh,u exp(u^2/2), etc...)
     -> (Double -> Double)          -- ^ derivative of transfer function
     -> NormType                    -- ^ type of normalisation: Infinity, PNorm1, PNorm2
     -> Double                      -- ^ convergence tolerance for feature vectors
     -> Matrix Double               -- ^ weight matrix
     -> [Matrix Double]             -- ^ input data in chunks
     -> Matrix Double               -- ^ ica transform (weight matrix)
ica' _ _  _ _ _ []     = error "no sample data"
ica' g g' n t w (x:xs) = let w' = normalise n $ decorrelate $ update g g' w x
                             in if converged t w w' 
                                then w'
                                else ica' g g' n t w' (xs ++ [x])

ica :: Int                         -- ^ random seed
    -> (Double -> Double)          -- ^ transfer function (tanh,u exp(u^2/2), etc...)
    -> (Double -> Double)          -- ^ derivative of transfer function
    -> NormType                    -- ^ type of normalisation: Infinity, PNorm1, PNorm2
    -> Double                      -- ^ convergence tolerance for feature vectors
    -> Int                         -- ^ output dimensions
    -> Int                         -- ^ sampling size (must be smaller than length of data)
    -> I.Array Int (Vector Double) -- ^ data
    -> (I.Array Int (Vector Double),Matrix Double) -- ^ transformed data, ica transform
ica r g g' n t o s a = let i = I.rangeSize $ I.bounds a
                           w = random_vector r (o,i)
                           x' = fromRows $ I.elems a
                           -- next line is BAD if distribution not stationary
                           x = concat $ toBlocksEvery i s x'
                           w' = ica' g g' n t w x
                           y = w' <> x'
                       in (I.listArray (1,o) $ toRows y,w') 

-----------------------------------------------------------------------------

-- | ICA with default values: no dimension reduction, euclidean norms, 16 sample groups, sigmoid
icaDefaults :: Int                         -- ^ random seed
            -> I.Array Int (Vector Double) -- ^ data
            -> (I.Array Int (Vector Double),Matrix Double) -- ^ transformed data, ica transform
icaDefaults r a = let c = I.rangeSize $ I.bounds a
                      s = (dim $ (a I.! 1)) `div` 16
                  in ica r sigmoid sigmoid' PNorm2 0.0000001 (c-1) s a

-----------------------------------------------------------------------------

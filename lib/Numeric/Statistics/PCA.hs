-----------------------------------------------------------------------------
-- |
-- Module      :  Numeric.Statistics.PCA
-- Copyright   :  (c) A. V. H. McPhail 2010, 2014, 2017
-- License     :  BSD3
--
-- Maintainer  :  haskell.vivian.mcphail <at> gmail <dot> com
-- Stability   :  provisional
-- Portability :  portable
--
-- Principal Components Analysis
--
-----------------------------------------------------------------------------

module Numeric.Statistics.PCA (
                               pca, pcaN, pcaTransform, pcaReduce, pcaReduceN
                          ) where


-----------------------------------------------------------------------------

import qualified Data.Array.IArray as I 

import Prelude hiding ((<>))
import Data.List(sortBy)
import Data.Ord(comparing)
import Numeric.LinearAlgebra
import Numeric.GSL.Statistics

--import Numeric.Statistics

-----------------------------------------------------------------------------

-- | find the principal components of multidimensional data greater than
--    the threshhold
pca :: I.Array Int (Vector Double)     -- the data
    -> Double                         -- eigenvalue threshold
    -> (Vector Double, Matrix Double) -- Eignevalues, Principal components
pca d q = let d' = fmap (\x -> x - (scalar $ mean x)) d -- remove the mean from each dimension
              d'' = fromColumns $ I.elems d'
              (_,vec',uni') = svd d''
              vec = toList vec'
              uni = toColumns uni'
              v' = zip vec uni
              v = filter (\(x,_) -> x > q) v'  -- keep only eigens > than parameter
              (eigs,vs) = unzip v
          in (fromList eigs,fromColumns vs) 

-- | find N greatest principal components of multidimensional data
--    according to size of the eigenvalue
pcaN :: I.Array Int (Vector Double)    -- the data
     -> Int                            -- number of components to return
     -> (Vector Double, Matrix Double) -- Eignevalues, Principal components
pcaN d n = let d' = fmap (\x -> x - (scalar $ mean x)) d -- remove the mean from each dimension
               d'' = fromColumns $ I.elems d'
               (_,vec',uni') = svd d''
               vec = toList vec'
               uni = toColumns uni'
               v' = zip vec uni
               v = take n $ reverse $ sortBy (comparing fst) v'
               (eigs,vs) = unzip v
           in (fromList eigs,fromColumns vs) 

v1 = fromList [1,2,3,4,5,6::Double]
v2 = fromList [2,3,4,5,6,7::Double]
v3 = fromList [3,4,5,6,7,8::Double]

a = fromColumns [v1,v2,v3]
b = I.listArray (1,3::Int) [v1,v2,v3] :: I.Array Int (Vector Double)
                
-- | perform a PCA transform of the original data (remove mean)
-- |     Final = M^T Data^T
pcaTransform :: I.Array Int (Vector Double)    -- ^ the data
             -> Matrix Double                  -- ^ the principal components
             -> I.Array Int (Vector Double)    -- ^ the transformed data
pcaTransform d m = let d' = fmap (\x -> x - (scalar $ mean x)) d -- remove the mean from each dimension
                   in I.listArray (1,cols m) $ toRows $ (tr' m) <> (fromRows $ I.elems d')

-- | perform a dimension-reducing PCA modification, 
--     using an eigenvalue threshhold
pcaReduce :: I.Array Int (Vector Double)      -- ^ the data
          -> Double                           -- ^ eigenvalue threshold
          -> I.Array Int (Vector Double)      -- ^ the reduced data
pcaReduce d q = let u = fmap (scalar . mean) d
                    d' = zipWith (-) (I.elems d) (I.elems u)
                    d'' = fromColumns d'
                    (_,vec',uni') = svd d''
                    vec = toList vec'
                    uni = toColumns uni'
                    v' = zip vec uni
                    v = filter (\(x,_) -> x > q) v'  -- keep only eigens > than parameter
                    m = fromColumns $ snd $ unzip v
                 in I.listArray (I.bounds d) $ zipWith (+) (toRows $ m <> (tr' m) <> fromRows d') (I.elems u) 

-- | perform a dimension-reducing PCA modification, using N components
pcaReduceN :: I.Array Int (Vector Double)      -- ^ the data
           -> Int                              -- ^ N, the number of components
           -> I.Array Int (Vector Double)      -- ^ the reduced data, with n principal components
pcaReduceN d n = let u = fmap (scalar . mean) d
                     d' = zipWith (-) (I.elems d) (I.elems u)
                     d'' = fromColumns d'
                     (_,vec',uni') = svd d''
                     vec = toList vec'
                     uni = toColumns uni'
                     v' = zip vec uni
                     v = take n $ reverse $ sortBy (comparing fst) v'
                     m = fromColumns $ snd $ unzip v
                  in I.listArray (I.bounds d) $ zipWith (+) (toRows $ m <> (tr' m) <> fromRows d') (I.elems u) 

-----------------------------------------------------------------------------

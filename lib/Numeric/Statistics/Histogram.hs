-----------------------------------------------------------------------------
-- |
-- Module      :  Numeric.Statistics.Histogram
-- Copyright   :  (c) A. V. H. McPhail 2010
-- License     :  GPL-style
--
-- Maintainer  :  haskell.vivian.mcphail <at> gmail <dot> com
-- Stability   :  provisional
-- Portability :  portable
--
-- create histograms from density functions
--
-----------------------------------------------------------------------------

module Numeric.Statistics.Histogram (
                                     cumulativeToHistogram
                                    , gaussianHistogram
                                  ) where


-----------------------------------------------------------------------------

import Data.Packed.Vector

import qualified Numeric.GSL.Histogram as H
--import qualified Numeric.GSL.Histogram2D as H2

import qualified Numeric.GSL.Distribution.Continuous as C

--import Numeric.LinearAlgebra.Algorithms
--import Numeric.LinearAlgebra.Interface()

-----------------------------------------------------------------------------

vectorToTuples = toTuples . toList
    where toTuples []         = error "need a minimum of two elements"
          toTuples [_]        = error "need a minimum of two elements"
          toTuples [x1,x2]    = [(x1,x2)]
          toTuples (x1:x2:xs) = (x1,x2) : (toTuples (x2:xs))

-----------------------------------------------------------------------------

cumulativeToHistogram :: (Double -> Double)  -- ^ the cumulative distribution function D(x <= X)
                   -> Vector Double          -- ^ the bins
                   -> H.Histogram            -- ^ the resulting histogram
cumulativeToHistogram f v = H.addListWeighted (H.emptyRanges v) $ map (\(x1,x2) -> ((x1 + x2) / 2.0,f x2 - f x1)) (vectorToTuples v)

gaussianHistogram :: Double              -- ^ mean
                  -> Double              -- ^ standard deviation
                  -> Vector Double       -- ^ the bins
                  -> H.Histogram         -- ^ the resulting histogram
gaussianHistogram u s = cumulativeToHistogram (\x -> C.density_1p C.Gaussian C.Lower s (x-u))

-----------------------------------------------------------------------------

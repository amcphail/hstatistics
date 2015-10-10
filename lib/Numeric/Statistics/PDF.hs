{-# LANGUAGE MultiParamTypeClasses,
             FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Numeric.Statistics.PDF
-- Copyright   :  (c) A. V. H. McPhail 2010, 2014
-- License     :  BSD3
--
-- Maintainer  :  haskell.vivian.mcphail <at> gmail <dot> com
-- Stability   :  provisional
-- Portability :  portable
--
-- Probability Distribution Function interface
--
-----------------------------------------------------------------------------

module Numeric.Statistics.PDF (
                               -- * PDF
                               PDFFunction,
                               PDF(..),
                               -- * Creation
                               pdfFromFunction
                          ) where


-----------------------------------------------------------------------------

--import qualified Data.Packed.Vector as V
import qualified Data.Vector.Storable as V

import qualified Numeric.GSL.Histogram as H
import qualified Numeric.GSL.Histogram2D as H2

import Foreign.Storable
--import Numeric.Statistics

-----------------------------------------------------------------------------

-- ^ a probability distribution function
data PDFFunction a = P_Func (a -> Double)     -- p(x)

-----------------------------------------------------------------------------

-- ^ a PDF interface
class PDF b a where
    -- ^ calculate a probability
    probability :: b -> V.Vector a -> V.Vector Double

instance Storable b => PDF (PDFFunction b) b where
    probability (P_Func f) = V.map f

instance PDF H.Histogram Double where
    probability = H.prob

instance PDF H2.Histogram2D (Double,Double) where
    probability = H2.probPaired

-----------------------------------------------------------------------------

-- | create a PDF from an arbtrary function f :-> [0,1]
pdfFromFunction :: (a -> Double) -> PDFFunction a
pdfFromFunction = P_Func


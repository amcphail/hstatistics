-----------------------------------------------------------------------------
-- |
-- Module      :  Numeric.Statistics.Viterbi.Discrete.Lazy
-- Copyright   :  (c) A. V. H. McPhail 2011
-- License     :  BSD3
--
-- Maintainer  :  haskell.vivian.mcphail <at> gmail <dot> com
-- Stability   :  provisional
-- Portability :  portable
--
-- Viterbi Algorithm
--
-----------------------------------------------------------------------------

module Numeric.Statistics.Viterbi.Discrete.Lazy (
                          ) where


-----------------------------------------------------------------------------

import qualified Data.Array.IArray as I 

import Numeric.LinearAlgebra

import Numeric.GSL.Statistics

import Numeric.Statistics

import System.Random

-----------------------------------------------------------------------------

newtype Viterbi a = Viterbi {
      _memory :: Int
    , _states :: Int
    




train :: Int       -- ^ sequence length
#pragma once

#include "rts/rtsmaths.h"
#include "rts/rtsmodel.h"
#include "rts/rtsregionmodel.h"

typedef rts::rtsModelBits<rts::ar1Covariance, glmmr::LinearPredictor> BitsAR;
typedef rts::rtsModelBits<rts::nngpCovariance, glmmr::LinearPredictor> BitsNNGP;
typedef rts::rtsModelBits<rts::ar1Covariance, rts::regionLinearPredictor> BitsARRegion;
typedef rts::rtsModelBits<rts::nngpCovariance, rts::regionLinearPredictor> BitsNNGPRegion;
typedef rts::rtsModel<rts::rtsModelBits<rts::ar1Covariance, glmmr::LinearPredictor> > ModelAR;
typedef rts::rtsModel<rts::rtsModelBits<rts::nngpCovariance, glmmr::LinearPredictor> > ModelNNGP;
typedef rts::rtsRegionModel<rts::rtsModelBits<rts::ar1Covariance, rts::regionLinearPredictor> > ModelARRegion;
typedef rts::rtsRegionModel<rts::rtsModelBits<rts::nngpCovariance, rts::regionLinearPredictor> > ModelNNGPRegion;

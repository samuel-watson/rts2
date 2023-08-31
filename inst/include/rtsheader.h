#pragma once

#include "rts/rtsmaths.h"
#include "rts/rtsmodel.h"
#include "rts/rtsregionmodel.h"

typedef rts::rtsModel<rts::rtsModelBits<rts::ar1Covariance, glmmr::LinearPredictor> > ModelAR;
typedef rts::rtsModel<rts::rtsModelBits<rts::nngpCovariance, glmmr::LinearPredictor> > ModelNNGP;
typedef rts::rtsRegionModel<rts::rtsModelBits<rts::ar1Covariance, rts::regionLinearPredictor> > ModelARRegionG;
typedef rts::rtsRegionModel<rts::rtsModelBits<rts::nngpCovariance, rts::regionLinearPredictor> > ModelNNGPRegionG;
typedef rts::rtsRegionModel<rts::rtsModelBits<rts::ar1Covariance, glmmr::LinearPredictor> > ModelARRegion;
typedef rts::rtsRegionModel<rts::rtsModelBits<rts::nngpCovariance, glmmr::LinearPredictor> > ModelNNGPRegion;

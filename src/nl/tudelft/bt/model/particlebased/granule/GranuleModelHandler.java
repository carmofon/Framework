/*
 * Created on Sep 15, 2004
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Generation - Code and Comments
 */
package nl.tudelft.bt.model.particlebased.granule;

import nl.tudelft.bt.model.apps.ModelHandler;
import nl.tudelft.bt.model.exceptions.InvalidValueException;
import nl.tudelft.bt.model.exceptions.ModelException;
import nl.tudelft.bt.model.exceptions.MultigridSystemNotSetException;
import nl.tudelft.bt.model.exceptions.SystemEditViolationException;
import nl.tudelft.bt.model.multigrid.MultigridVariable;
import nl.tudelft.bt.model.multigrid.boundary_conditions.GranuleBoundaryConditions;
import nl.tudelft.bt.model.multigrid.boundary_layers.SphericalDilationBoundaryLayer;

/**
 * @author jxavier
 * 
 * TODO To change the template for this generated type comment go to Window -
 * Preferences - Java - Code Generation - Code and Comments
 */
public abstract class GranuleModelHandler extends ModelHandler {

	/**
	 * 
	 */
	public GranuleModelHandler() {
		super();
	}

	protected void builBiomassParticleContainer() {
		_m.buildGranuleBiomassContainer(_relativeMaximumRadius * _systemSize,
				_relativeMinimumRadius * _systemSize, _kShoving, FSHOVING);
	}


	protected void inoculateRandomly(int[] nCells) {
		if (!_physicalSystemIsEditable) {
			throw (new SystemEditViolationException(
					"Tried to innoculate system" + " when not in edit mode"));
		}
		_m.inoculateRandomlyMultispeciesInGranule(nCells);
	}
	
	protected void inoculateRandomlyInGranuleAtRadius (int[] nCells, float r) {
		if (!_physicalSystemIsEditable) {
			throw (new SystemEditViolationException(
					"Tried to innoculate system" + " when not in edit mode"));
		}
		_m.inoculateRandomlyMultispeciesInGranuleAtRadius(nCells, r);
	}

	protected void inoculateRandomlyInsideRadius(int[] nCells, float r) {
		if (!_physicalSystemIsEditable) {
			throw (new SystemEditViolationException(
					"Tried to innoculate system" + " when not in edit mode"));
		}
		_m.inoculateRandomlyMultispeciesInColony(nCells, r);
	}
	
	protected void inoculateRandomlyAtRadius(int[] nCells, float r) {
		if (!_physicalSystemIsEditable) {
			throw (new SystemEditViolationException(
					"Tried to innoculate system" + " when not in edit mode"));
		}
		_m.inoculateRandomlyMultispeciesAtRadius(nCells, r);
	}

	protected void inoculateRandomlyEverywhere(int[] nCells) {
		if (!_physicalSystemIsEditable) {
			throw (new SystemEditViolationException(
					"Tried to innoculate system" + " when not in edit mode"));
		}
		_m.inoculateRandomlyMultispeciesEverywhere(nCells);
	}

	
	
	
	@Override
	public void initializeDiffusionReactionSystem() throws ModelException {
		super.initializeDiffusionReactionSystem();
		// This overrides the max height so that we can implement cyclic boundaries
		_m.setVerticalCutoffSize(Float.POSITIVE_INFINITY);
	}

	@Override
	protected void inoculate() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void initializeDetachmentFunction() {
		// TODO Auto-generated method stub
		
	}

	/**
	 * Turns off cell shedding
	 */
	protected void turnShedingOn() {
		GranuleBiomassParticle.turnSheddingOn();
	}
}
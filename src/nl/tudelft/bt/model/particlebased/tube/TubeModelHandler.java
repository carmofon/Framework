/*
 * Created on Sep 15, 2004
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Generation - Code and Comments
 */
package nl.tudelft.bt.model.particlebased.tube;

import java.util.Collection;
import java.util.Iterator;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.apps.ModelHandler;
import nl.tudelft.bt.model.exceptions.ModelException;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.exceptions.MultigridSystemNotSetException;
import nl.tudelft.bt.model.exceptions.SystemEditViolationException;
import nl.tudelft.bt.model.multigrid.MultigridVariable;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;
import nl.tudelft.bt.model.multigrid.TubeVelocityField;
import nl.tudelft.bt.model.multigrid.boundary_conditions.TubeBoundaryConditions;

/**
 * @author jxavier
 * 
 * TODO To change the template for this generated type comment go to Window -
 * Preferences - Java - Code Generation - Code and Comments
 */
public abstract class TubeModelHandler extends ModelHandler {
	private float _tubeRadius;

	private float _tubeLength;
	
	private float _volumetricFlowRate;

	private float _viscosityLiquid;

	private float _referenceDiffusivity;
	
	/**
	 * Crerate a new instance of TubeModelHandler with the radius
	 * @param tubeRadius
	 * @param volumetricFlowRate
	 * @param viscosity
	 *            the value of the liquid's viscosity
	 */
	public TubeModelHandler(float tubeRadius, float tubeLength, float volumetricFlowRate,
			float viscosity) {
		super();
		_tubeRadius = tubeRadius;
		_tubeLength = tubeLength;
		_volumetricFlowRate = volumetricFlowRate;
		_viscosityLiquid = viscosity;
	}

	protected void builBiomassParticleContainer() {
		// check if the system is 2D
		if (_geometry != 2) {
			throw new ModelRuntimeException("Tube reactors may only be 2D");
		}
		// create the biomass container with tube geometry
		_m.buildTubeBiomassContainer(_relativeMaximumRadius * _systemSize,
				_relativeMinimumRadius * _systemSize, _kShoving, FSHOVING,
				_tubeRadius);
	}

	protected void createBoundaryLayer(float h)
			throws MultigridSystemNotSetException {
		_boundaryLayer = TubeVelocityField.getInstance().getBoundaryLayer();
		//create new boudary conditions
		MultigridVariable.setBoundaryConditions(new TubeBoundaryConditions(
				_tubeRadius));
	}

	protected void inoculateRandomly(int[] nCells) {
		if (!_physicalSystemIsEditable) {
			throw (new SystemEditViolationException(
					"Tried to innoculate system" + " when not in edit mode"));
		}
		_m.inoculateRandomlyMultispeciesInTube(nCells);
	}

	/* (non-Javadoc)
	 * @see nl.tudelft.bt.model.apps.ModelHandler#initializeSoluteConcentrations()
	 */
	public void initializeSoluteConcentrations()
			throws MultigridSystemNotSetException {
		computeFluidVelocityInTube();
		super.initializeSoluteConcentrations();
	}
	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.apps.ApplicationComponent#performFullIteration()
	 */
	public void detach() throws ModelException {
		computeFluidVelocityInTube();
		super.detach();
	}

	private void computeFluidVelocityInTube()
			throws MultigridSystemNotSetException {
		//get the reference diffusivity (average of all solutes)
		Collection sp = Model.model().getSoluteSpecies();
		float sumDiffusivities = 0;;
		for (Iterator iter = sp.iterator(); iter.hasNext();) {
			SoluteSpecies element = (SoluteSpecies) iter.next();
			sumDiffusivities += element.getDiffusivity();
		}
		_referenceDiffusivity = sumDiffusivities/sp.size();
		//
		TubeVelocityField.getInstance().setFlowProperties(_volumetricFlowRate,
				_viscosityLiquid, _tubeRadius, _tubeLength,
				_referenceDiffusivity);
		TubeVelocityField.getInstance().solveVelocityFieldAndComputeShear();
	}
}
package nl.tudelft.bt.model.work.quorumsensing;

import nl.tudelft.bt.model.multigrid.MultigridVariable;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;
import nl.tudelft.bt.model.reaction.ProcessFactor;

/**
 * Implements a step function for simulation of quorum sensing switch.
 * Turns rate off once concentration of autoinducer is above QSthreshold
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu)
 * @see UpRegulationStep
 */
public class DownRegulationStep extends ProcessFactor {
	protected MultigridVariable _autoinducer; 
	private float _QSthreshold;

	/**
	 * @param 
	 */
	public DownRegulationStep(MultigridVariable autoinducer, float QSthreshold) {
		_autoinducer = autoinducer;
		_QSthreshold = QSthreshold;
		
	}

	public float getValue() {
		float conc = _autoinducer.getValue();
		conc = (conc < 0 ? 0 : conc);
		if (conc > _QSthreshold) {
			return 0;
		}
		return 1;
	}

	public float getMaximumValue() {
		return 1;
	}

	public float getDerivative(SoluteSpecies c) {
		return 0f;
	}

}

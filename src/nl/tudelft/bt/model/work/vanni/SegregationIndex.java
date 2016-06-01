package nl.tudelft.bt.model.work.vanni;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.apps.output.VariableSeries;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.reaction.ProcessFactor;

import java.util.Collection;
import java.util.Iterator;

public class SegregationIndex extends VariableSeries {
	private BacteriocinBiomassSpecies _sensitive;
	private BacteriocinBiomassSpecies _producer;
	protected ProcessFactor _mO;
	float _Lbac;
	protected MutatorBiomassSpecies _species;
	
	public SegregationIndex(BacteriocinBiomassSpecies sensitive,
			BacteriocinBiomassSpecies producer, float Lbac) {
		super("SegragationIndex", "Time [h]",
				"SI[adimensional]");
		setX(Model.model().getTimeSeries());
		_producer = producer;
		_sensitive = sensitive;
		_mO = producer._m0;
		_Lbac = Lbac;
	}
	
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.model.apps.output.VaribleSeries#getLastY()
	 */
	public float getLastY() {
		float v = 0;
		float Sum_mO =0;
		Collection list = Model.model().biomassContainer
				.getBiomassAsBiomassParticleCollection();
		for (Iterator iter = list.iterator(); iter.hasNext();) {
			// get the current bacterium
			BiomassParticle focal = (BiomassParticle) iter.next();
			if (focal.isOfSpecies(_producer)) {
				for (Iterator iter2 = list.iterator(); iter2.hasNext();) {
					BiomassParticle other = (BiomassParticle) iter2.next();
					if (other.isOfSpecies(_sensitive)) {
						v += _Lbac * _mO.getValue() * 1 / focal.distanceTo(other);
					}
				}
			}
		}
		// NEED TO AVERAGE IT AMONG ALL THE ACTIVELY GROWING PRODUCERS CELLS
		for (Iterator iter = list.iterator(); iter.hasNext();){
				// get the current bacterium
				BiomassParticle focal = (BiomassParticle) iter.next();
				if (focal.isOfSpecies(_producer)) {
					Sum_mO += _mO.getValue(); 
		}
		}
		//
		int sizeX = getXArray().getSize();
		int sizeY = getYArray().getSize();
		if (sizeY < sizeX) {
			getYArray().add(v/Sum_mO);
		}
		// every time getY is invoked, the array is updated
		return super.getLastY();
	}

}


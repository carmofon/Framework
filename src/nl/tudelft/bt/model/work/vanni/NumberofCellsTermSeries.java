package nl.tudelft.bt.model.work.vanni;

import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.apps.output.VariableSeries;
import nl.tudelft.bt.model.particlebased.BiomassParticle;

public class NumberofCellsTermSeries extends VariableSeries {
	private BiomassSpecies _cheater;
	private BiomassSpecies _producer;

	public NumberofCellsTermSeries(BiomassSpecies cheater,
			BiomassSpecies producer) {
		super("Relatedness measure", "Time [h]",
				"Ncheater/Nproducers");
		setX(Model.model().getTimeSeries());
		_cheater = cheater;
		_producer = producer;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.model.apps.output.VaribleSeries#getLastY()
	 */
	public float getLastY() {
		float v1 = 0;
		float v2 = 0;
		Collection list = Model.model().biomassContainer
				.getBiomassAsBiomassParticleCollection();
		for (Iterator iter = list.iterator(); iter.hasNext();) {
			// get the current bacterium
			BiomassParticle focal = (BiomassParticle) iter.next();
			//Calculate term 1
			if (focal.isOfSpecies(_producer)){
				v1 += 1;
			}else{
				v2 += 1;
			}				
		}
		//
		int sizeX = getXArray().getSize();
		int sizeY = getYArray().getSize();
		if (sizeY < sizeX) {
			getYArray().add(v2/v1);
		}
		// every time getY is invoked, the array is updated
		return super.getLastY();
	}

}

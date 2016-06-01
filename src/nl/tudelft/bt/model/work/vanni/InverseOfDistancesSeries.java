package nl.tudelft.bt.model.work.vanni;

import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.apps.output.VariableSeries;
import nl.tudelft.bt.model.particlebased.BiomassParticle;

public class InverseOfDistancesSeries extends VariableSeries {
	private BiomassSpecies _sensitive;
	private BiomassSpecies _producer;

	public InverseOfDistancesSeries(BiomassSpecies sensitive,
			BiomassSpecies producer) {
		super("Sum of inverse distances", "Time [h]",
				"Sum of inverse distances [1/um]");
		setX(Model.model().getTimeSeries());
		_producer = producer;
		_sensitive = sensitive;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.model.apps.output.VaribleSeries#getLastY()
	 */
	public float getLastY() {
		float v = 0;
		Collection list = Model.model().biomassContainer
				.getBiomassAsBiomassParticleCollection();
		for (Iterator iter = list.iterator(); iter.hasNext();) {
			// get the current bacterium
			BiomassParticle focal = (BiomassParticle) iter.next();
			if (focal.isOfSpecies(_sensitive)) {
				for (Iterator iter2 = list.iterator(); iter2.hasNext();) {
					BiomassParticle other = (BiomassParticle) iter2.next();
					if (other.isOfSpecies(_producer)) {
						v += 1 / focal.distanceTo(other);
					}
				}
			}
		}
		//
		int sizeX = getXArray().getSize();
		int sizeY = getYArray().getSize();
		if (sizeY < sizeX) {
			getYArray().add(v);
		}
		// every time getY is invoked, the array is updated
		return super.getLastY();
	}

}

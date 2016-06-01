package nl.tudelft.bt.model.work.vanni;

import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.apps.output.VariableSeries;
import nl.tudelft.bt.model.particlebased.BiomassParticle;

public class ComputeRelatednessTermSeries extends VariableSeries {
	private BiomassSpecies _cheater;
	private BiomassSpecies _producer;

	public ComputeRelatednessTermSeries(BiomassSpecies cheater,
			BiomassSpecies producer) {
		super("Relatedness measure", "Time [h]",
				"Volume per Inverse Distance [um^3/um]");
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
		float v3 = 0;
		Collection list = Model.model().biomassContainer
				.getBiomassAsBiomassParticleCollection();
		for (Iterator iter = list.iterator(); iter.hasNext();) {
			// get the current bacterium
			BiomassParticle focal = (BiomassParticle) iter.next();
	
			//Calculate term 1
			if (focal.isOfSpecies(_producer)){
				for (Iterator iter2 = list.iterator(); iter2.hasNext();) {
					BiomassParticle other = (BiomassParticle) iter2.next();
					if (other.isOfSpecies(_producer)) {
						if (other.equals(focal) == false){	
							float tmp= focal.distanceTo(other);
							if (Model.model().getDimensionality() == 3) {
								v1 += 4* Math.PI * Math.pow(other.getRadius(), 3)/tmp;
							}
							else {
								// calculate the cylinder volume
								v1 += 2* Math.PI *Math.pow(other.getRadius(), 2)/tmp;
							}
						}
					}
				}
			}
				
			//Calculate term 2
			if (focal.isOfSpecies(_producer)){
				if (Model.model().getDimensionality() == 3) {
					v2 += 4* Math.PI * Math.pow(focal.getRadius(), 3)/focal.getRadius();
				}
				else {
					// calculate the cylinder volume
					v2 += 2* Math.PI * Math.pow(focal.getRadius(), 2)/focal.getRadius();
				}
			}
			
			//Calculate term 3
			if (focal.isOfSpecies(_cheater)){
				for (Iterator iter2 = list.iterator(); iter2.hasNext();) {
					BiomassParticle other = (BiomassParticle) iter2.next();
					if (other.isOfSpecies(_producer)) {
						if (Model.model().getDimensionality() == 3) {
							v3 += 4* Math.PI * Math.pow(other.getRadius(), 3)/focal.distanceTo(other);
						}
						else {
							// calculate the cylinder volume
							v3 += 2* Math.PI * Math.pow(other.getRadius(), 2)/ focal.distanceTo(other);
						}
					}
				}
			}
		}
		//
		int sizeX = getXArray().getSize();
		int sizeY = getYArray().getSize();
		if (sizeY < sizeX) {
			getYArray().add(v1+v2-v3);
		}
		// every time getY is invoked, the array is updated
		return super.getLastY();
	}

}

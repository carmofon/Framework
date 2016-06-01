package nl.tudelft.bt.model.work.quorumsensing;

import java.awt.Color;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.multigrid.MultigridVariable;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;

/**
 * Implements the mutator particle
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Oct 18, 2006
 */
public class QSBiomassParticle extends BiomassParticle {
	public QSBiomassParticle(QSBiomassSpecies s) {
		super(s);
	}

	//TODO check if these methods can be removed
	@Override
	public BiomassParticle divide() {
		QSBiomassParticle daughter = (QSBiomassParticle) super
				.divide();
		return daughter;
	}

	@Override
	protected Object clone() throws CloneNotSupportedException {
		return super.clone();
	}

	@Override
	public Color getColorCore() {
		//The problem was that setCurrentLocation was calle after getValue!
		MultigridVariable.setCurrentLocation(getCenter());
		double a = ((QSBiomassSpecies)_biomassSpecies)._autoInducer.getValue();
		double t = ((QSBiomassSpecies)_biomassSpecies)._threshold;
		Color cBelow = ((QSBiomassSpecies)_biomassSpecies)._colorBelowQuorum;
		Color cAbove = ((QSBiomassSpecies)_biomassSpecies)._colorAboveQuorum;
		
		// TODO Auto-generated method stub
		return (a > t ? cAbove : cBelow);
	}
}

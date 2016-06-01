package nl.tudelft.bt.model.work.quorumsensing;

import java.awt.Color;

import nl.tudelft.bt.model.ContinuousCoordinate;
import nl.tudelft.bt.model.exceptions.MultigridSystemNotSetException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;

public class QuorumSensingSpecies extends ParticulateSpecies {
	private double _threshold;
	private SoluteSpecies _autoInducer;
	private boolean _quorum;
	private Color _cBelowQuorum;
	private Color _cAboveQuorum;
	
	public QuorumSensingSpecies(String n, float dens, SoluteSpecies autoInducer, double threshold, 
			Color cBelowQuorum, Color cAboveQuorum) throws MultigridSystemNotSetException {
		super(n, dens, cBelowQuorum);
		_cBelowQuorum = cBelowQuorum;
		_cAboveQuorum = cAboveQuorum;
		_autoInducer = autoInducer;
		_threshold = threshold;
		
	}

	@Override
	public float getSpecificRate(ContinuousCoordinate c) {
		float r = super.getSpecificRate(c);
		float a = _autoInducer.getValue();
		if (a > _threshold)
			quorumOn();
		else
			quorumOff();
		return r;
	}

	@Override
	public Color getColor() {
		return (_quorum ? _cAboveQuorum : _cBelowQuorum);
	}
	
	private void quorumOff() {
		_quorum = false;
	}
	
	private void quorumOn () {
		_quorum = true;
	}
}

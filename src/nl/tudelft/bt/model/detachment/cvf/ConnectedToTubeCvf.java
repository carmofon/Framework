/*
 * Created on Sep 14, 2004
 */
package nl.tudelft.bt.model.detachment.cvf;

import nl.tudelft.bt.model.particlebased.tube.TubeBiomassParticleContainer;

/**
 * @author jxavier
 */
public class ConnectedToTubeCvf extends ConnectedToTopCvf {
	private TubeBiomassParticleContainer _container;

	/**
	 * @param n
	 * @param m
	 * @param reference
	 *            to the tube biomass container
	 */
	public ConnectedToTubeCvf(int n, int m, TubeBiomassParticleContainer c) {
		//note: tube reactor only defined for 2D
		super(n, m, 1);
		_container = c;
	}

	/**
	 * The initiator starts with the top values
	 */
	protected void initializeCvf() {
		//initialize center to false
		for (int i = 0; i < _n; i++) {
			for (int j = 0; j < _m; j++) {
				for (int k = 0; k < _l; k++) {
					if (_container.belongsToTubeBorder(i, j, k)) {
						_cvf[i][j][k] = _matrixToFilter[i][j][k];
					} else
						_cvf[i][j][k] = false;
				}
			}
		}
	}
}
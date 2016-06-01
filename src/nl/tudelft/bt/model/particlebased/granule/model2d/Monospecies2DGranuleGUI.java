/*
 * Created on Sep 8, 2004
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Generation - Code and Comments
 */
package nl.tudelft.bt.model.particlebased.granule.model2d;

import java.awt.*;
import java.awt.event.*;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Hashtable;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import nl.tudelft.bt.model.exceptions.GrowthProblemException;
import nl.tudelft.bt.model.exceptions.ModelException;

/**
 * @author jxavier
 * 
 * TODO To change the template for this generated type comment go to Window -
 * Preferences - Java - Code Generation - Code and Comments
 */
public class Monospecies2DGranuleGUI extends JPanel implements ActionListener {
	private JButton _buttonRun, _buttonPause, _buttonRestart;

	JPanel _buttonPane;

	private boolean _pauseSwitch;

	private static final String RUN = "  run   ";

	private static final String PAUSE = " pause  ";

	private static final String CONTINUE = "continue";

	private static final String RESTART = "restart ";

	private Monospecies2DGranuleModel _m;

	private Runner _runner;

	private ContourPlotHolder cp;

	private SliderCreator _uMaxRangeSlider;

	private SliderCreator _bDecayRangeSlider;

	private SliderCreator _inocolumRangeSlider;

	private SliderCreator _cORangeSlider;

	private SliderCreator _boundaryLayerRangeSlider;

	private SliderCreator _detachmentRangeSlider;

	/**
	 * Thread that runs the simulation
	 * 
	 * @author jxavier
	 */
	private class Runner extends Thread {
		public void run() {
			try {
				_m.start();
			} catch (Exception e) {
				if (e instanceof GrowthProblemException) {
					// biofilm is gone, so return
					return;
				}
				e.printStackTrace();
			}
		}
	}

	private class SliderCreator implements ChangeListener {
		float _max;

		float _min;

		float _init;

		String _label;

		String _units;

		JSlider _slider;

		private JLabel _sliderLabel;

		private JPanel _comp;

		private NumberFormat _frmt;

		/**
		 * Create a new slider
		 * 
		 * @param label
		 *            label to display
		 * @param units
		 *            units for the value read by the slider
		 * @param min
		 *            minimum in the scale
		 * @param max
		 *            maximum in the scale
		 * @param init
		 *            initial value
		 * @param frmt
		 *            the number format
		 */
		public SliderCreator(String label, String units, float min, float max,
				float init, NumberFormat frmt) {
			_label = label;
			_units = units;
			_max = max;
			_min = min;
			_init = init;
			_frmt = frmt;
			// Create the slider.
			_slider = new JSlider(JSlider.HORIZONTAL, 0, 100, 50);
			_slider.addChangeListener(this);
			_slider.setPaintLabels(true);
			// Create the label.
			_sliderLabel = new JLabel(label, JLabel.CENTER);
			setSliderLabel(getValue());
			Font f = _sliderLabel.getFont();
			_sliderLabel.setFont(new Font(f.getName(), f.getStyle(), f
					.getSize() - 3));
			_sliderLabel.setAlignmentX(Component.CENTER_ALIGNMENT);
		}

		JLabel createLabel(float v) {
			JLabel lab = new JLabel(_frmt.format(v));
			Font f = lab.getFont();
			lab.setFont(new Font(f.getName(), f.getStyle(), f.getSize() - 3));
			return lab;
		}

		/**
		 * Return the slider, label and text field organized in a container
		 * 
		 * @return
		 */
		public JComponent getComponent() {
			if (_comp == null) {
				// set the initial value
				_slider.setValue(invertValue(_init));
				// create the labels
				createSliderTickLabels();
				// create the basis component
				_comp = new JPanel();
				_comp.setLayout(new BoxLayout(_comp, BoxLayout.X_AXIS));
				// create the component for slider and label
				JPanel sliderAndLabelPanel = new JPanel();
				sliderAndLabelPanel.setLayout(new BoxLayout(
						sliderAndLabelPanel, BoxLayout.Y_AXIS));
				sliderAndLabelPanel.add(_sliderLabel);
				sliderAndLabelPanel.add(_slider);
				_comp.add(sliderAndLabelPanel);
				// add the text field
			}
			return _comp;
		}

		void createSliderTickLabels() {
			// Turn on labels at major tick marks.
			_slider.setMajorTickSpacing(50);
			_slider.setMinorTickSpacing(25);
			_slider.setPaintTicks(true);
			// Create the label table
			// NOTE: must be performed after construction!
			Hashtable labelTable = new Hashtable();
			labelTable.put(new Integer(0), createLabel(_min));
			labelTable.put(new Integer(25), createLabel(computeValue(25)));
			labelTable.put(new Integer(50), createLabel(computeValue(50)));
			labelTable.put(new Integer(75), createLabel(computeValue(75)));
			labelTable.put(new Integer(100), createLabel(_max));
			_slider.setLabelTable(labelTable);
		}

		/**
		 * @return the value slected
		 */
		public float getValue() {
			return computeValue(_slider.getModel().getValue());
		}

		/**
		 * Coompute the value for the slider
		 * 
		 * @param v
		 * @return
		 */
		float computeValue(int v) {
			return (float) (v) / 100f * (_max - _min) + _min;
		}

		int invertValue(float v) {
			return (int) (100f * (v - _min) / (_max - _min));
		}

		public void stateChanged(ChangeEvent arg0) {
			setSliderLabel(getValue());
		}

		private void setSliderLabel(float v) {
			_sliderLabel.setText("<html>" + _label + ": "
					+ _frmt.format(getValue()) + " " + _units + "</html>");
		}
	}

	private class LogSliderCreator extends SliderCreator {
		private float _gamma;

		/**
		 * @param label
		 * @param min
		 * @param max
		 * @param init
		 */
		public LogSliderCreator(String label, String units, float min,
				float max, float init, NumberFormat frmt, float gamma) {
			super(label, units, min, max, init, frmt);
			_gamma = gamma;
		}

		float computeValue(int v) {
			return computeGammaValue(v) * (_max - _min) + _min;
		}

		int invertValue(float v) {
			return (int) (Math.pow(((v - _min) / (_max - _min)), 1 / _gamma) * 100f);
		}

		private float computeGammaValue(int v) {
			return (float) Math.pow(((float) v) / 100f, _gamma);
		}

		void createSliderTickLabels() {
			int aux1, aux2;
			// Turn on labels at major tick marks.
			_slider.setMajorTickSpacing(50);
			_slider.setPaintTicks(true);
			// Create the label table
			// NOTE: must be performed after construction!
			Hashtable labelTable = new Hashtable();
			labelTable.put(new Integer(0), createLabel(_min));
			labelTable.put(new Integer(50), createLabel(computeValue(50)));
			labelTable.put(new Integer(100), createLabel(_max));
			_slider.setLabelTable(labelTable);
		}
	}

	/**
	 * @throws java.awt.HeadlessException
	 */
	public Monospecies2DGranuleGUI() {
		// panes will be set vertically
		setLayout(new BoxLayout(this, BoxLayout.X_AXIS));

		// The left pane
		JPanel leftPane = new JPanel();
		leftPane.setLayout(new BorderLayout());
		// creat the plot
		cp = new ContourPlotHolder();
		leftPane.add(cp.panel, BorderLayout.CENTER);

		// create a pane with the 3 buttons
		_buttonPane = new JPanel();

		_buttonRun = new JButton(RUN);
		_buttonRun.setActionCommand(RUN);
		_buttonRun.addActionListener(this);
		_buttonRun.setToolTipText("Start simulation");

		_buttonPause = new JButton(PAUSE);
		_buttonPause.setActionCommand(PAUSE);
		_buttonPause.addActionListener(this);
		_buttonPause.setToolTipText("Pause/continue simulation");

		_buttonRestart = new JButton(RESTART);
		_buttonRestart.setActionCommand(RESTART);
		_buttonRestart.addActionListener(this);
		_buttonRestart.setToolTipText("restart/re-inoculate system");

		_buttonPane.add(_buttonRun);
		_buttonPane.add(_buttonPause);
		_buttonPane.add(_buttonRestart);
		leftPane.add(_buttonPane, BorderLayout.SOUTH);

		// initialize the pause switch
		_buttonRun.setEnabled(true);
		_buttonPause.setEnabled(false);
		_buttonRestart.setEnabled(true);
		_pauseSwitch = false;

		// The right pane
		JPanel rightPane = new JPanel();
		rightPane.setLayout(new BoxLayout(rightPane, BoxLayout.Y_AXIS));
		// second pane - properties of the microorganism
		JPanel sliderPane1 = new JPanel();
		sliderPane1.setBorder(BorderFactory
				.createTitledBorder("Properties of microorganism"));
		sliderPane1.setLayout(new BoxLayout(sliderPane1, BoxLayout.Y_AXIS));
		// the qMax Slider
		_uMaxRangeSlider = new SliderCreator("u<sup>max</sup>",
				"h<sup>-1</sup>", 0.001f, 0.6f, 0.53f, new DecimalFormat(
						"0.###"));
		sliderPane1.add(_uMaxRangeSlider.getComponent());
		// the Ko Slider
		_bDecayRangeSlider = new SliderCreator("b<sub>decay</sub>",
				"h<sup>-1</sup>", 0, 0.01f, 0f, new DecimalFormat("0.#####E0"));
		sliderPane1.add(_bDecayRangeSlider.getComponent());
		// the inoculum Slider
		_inocolumRangeSlider = new SliderCreator("number of initial particles",
				"", 1, 500, 200, new DecimalFormat("0"));
		sliderPane1.add(_inocolumRangeSlider.getComponent());
		// add the second pane
		rightPane.add(sliderPane1);

		// second pane - properties of the microorganism
		JPanel sliderPane2 = new JPanel();
		sliderPane2.setBorder(BorderFactory
				.createTitledBorder("Environment properties"));
		sliderPane2.setLayout(new BoxLayout(sliderPane2, BoxLayout.Y_AXIS));
		// the Co Slider
		_cORangeSlider = new SliderCreator("C_O", "g<sub>O</sub>/L", 9e-4f,
				9e-3f, 9e-3f, new DecimalFormat("0.#####E0"));
		sliderPane2.add(_cORangeSlider.getComponent());
		// the Boundary layer Slider
		_boundaryLayerRangeSlider = new SliderCreator(
				"Boundary layer thickness", "x10<sup>-6</sup>m", 0, 400f, 100f,
				new DecimalFormat("0.###"));
		sliderPane2.add(_boundaryLayerRangeSlider.getComponent());
		// the Detachment layer Slider
		_detachmentRangeSlider = new LogSliderCreator("Detachment constant",
				"x10<sup>-15</sup> gCOD-H.um<sup>-4<sup>.h<sup>-1<sup>", 0.1f,
				0.7f, 0.6967f, new DecimalFormat("0.##E0"), 4f);
		sliderPane2.add(_detachmentRangeSlider.getComponent());
		// add the third pane
		rightPane.add(sliderPane2);

		// add the left and right panes
		add(leftPane);
		add(rightPane);

		try {
			// initialize model and reset
			_m = new Monospecies2DGranuleModel(this);
		} catch (ModelException e) {
			e.printStackTrace();
		}
		updateParticleVisualizer();
	}

	public void actionPerformed(ActionEvent e) {
		if (RUN.equals(e.getActionCommand())) {
			_buttonRun.setEnabled(false);
			_buttonRestart.setEnabled(true);
			_buttonPause.setEnabled(true);
			if (_runner != null)
				if (_runner.isAlive()) {
					// simulation is running, return
					return;
				}
			getRootPane().setDefaultButton(_buttonPause);
			_runner = new Runner();
			_runner.start();
		} else if (PAUSE.equals(e.getActionCommand())) {
			_pauseSwitch = !_pauseSwitch;
			if (_pauseSwitch) {
				_buttonPause.setText(CONTINUE);
			} else {
				_buttonPause.setText(PAUSE);
			}
			getRootPane().setDefaultButton(_buttonPause);
			// pause command to model
			_m.pause();
		} else if (RESTART.equals(e.getActionCommand())) {
			// restart command
			_pauseSwitch = false;
			_buttonPause.setText(PAUSE);
			_buttonPause.setEnabled(false);
			_buttonRestart.setEnabled(false);
			try {
				// Add command to reset model here
				_m.restart();
			} catch (ModelException e1) {
				e1.printStackTrace();
			}
			updateParticleVisualizer();
			_buttonRestart.setEnabled(true);
			_buttonRun.setEnabled(true);
			getRootPane().setDefaultButton(_buttonRun);
		}
	}

	public void stop() {
		// Stop the running thread.
		_runner = null;
	}

	public void stateChanged(ChangeEvent arg0) {
		// do nothing
	}

	public float getUserUMax() {
		return _uMaxRangeSlider.getValue();
	}

	public float getUserBDecay() {
		return _bDecayRangeSlider.getValue();
	}

	public int getUserInoculum() {
		return (int) _inocolumRangeSlider.getValue();
	}

	public float getUserCO() {
		return _cORangeSlider.getValue();
	}

	public float getUserBoundaryLayerThickness() {
		return _boundaryLayerRangeSlider.getValue();
	}

	public float getUserDetachment() {
		return _detachmentRangeSlider.getValue();
	}

	public static void main(String[] args) {
		// Create a frame and add the panel
		JFrame frm = new JFrame();
		Monospecies2DGranuleGUI gui = new Monospecies2DGranuleGUI();
		frm.getContentPane().add(gui);
		gui.getRootPane().setDefaultButton(gui._buttonRun);
		frm.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frm.pack();
		frm.setVisible(true);
	}

	public void updateParticleVisualizer() {
		cp.panel.updateParticleData();
	}
}
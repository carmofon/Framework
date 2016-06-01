package nl.tudelft.bt.model.gui;

import nl.tudelft.bt.model.apps.*;
import nl.tudelft.bt.model.exceptions.ModelException;

import javax.swing.*; //This is the final package name.
//import com.sun.java.swing.*; //Used by JDK 1.2 Beta 4 and all
//Swing releases before Swing 1.1 Beta 3.
import java.awt.*;
import java.awt.event.*;

/**
 * Window with bottons to control simultion (still in test phase)
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class Visualizer2DControler {
	private static ApplicationComponent _app;
	private JFrame _frm;

	public Component createComponents() {
		// button for requesting a full iteration
		JButton button1 = new JButton("Full iteration");
		button1.setMnemonic(KeyEvent.VK_I);
		button1.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					_app.performFullIteration();
				} catch (ModelException ex) {
				}
			}
		});

		// button for requesting a growth and division step
		JButton button2 = new JButton("grow and divide");
		button2.setMnemonic(KeyEvent.VK_I);
		button2.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					_app.performGrowthAndDivision();
				} catch (ModelException ex) {
					System.out.println(ex);
				}
			}
		});

		// button for requesting a full shoving relaxation
		JButton button3 = new JButton("spread completely");
		button3.setMnemonic(KeyEvent.VK_I);
		button3.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				_app.spreadCompletely();
			}
		});

		// button for requesting a full shoving relaxation
		JButton button4 = new JButton("shove completely");
		button4.setMnemonic(KeyEvent.VK_I);
		button4.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				_app.spreadByShovingCompletely();
			}
		});

		// button for requesting a full shoving relaxation
		JButton button5 = new JButton("spread by pressure 1");
		button5.setMnemonic(KeyEvent.VK_I);
		button5.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				_app.performSpreadingByPulling();
			}
		});

		// button for requesting a full shoving relaxation
		JButton button6 = new JButton("spread by pressure 2");
		button6.setMnemonic(KeyEvent.VK_I);
		button6.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				_app.performSpreadingByPulling();
			}
		});

		// button for requesting a shoving step
		JButton button7 = new JButton("single shoving step");
		button7.setMnemonic(KeyEvent.VK_I);
		button7.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				_app.performSpreadingStep();
			}
		});

		// button for requesting a shoving step
		JButton button8 = new JButton("continue iterating");
		button8.setMnemonic(KeyEvent.VK_I);
		button8.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				_app.continueIterating();
			}
		});

		// button for requesting a shoving step
		JButton button9 = new JButton("break iteration");
		button9.setMnemonic(KeyEvent.VK_I);
		button9.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				_app.stopIterating();
			}
		});

		/*
		 * An easy way to put space between a top-level container
		 * and its contents is to put the contents in a JPanel
		 * that has an "empty" border.
		 */
		JPanel pane = new JPanel();
		pane.setBorder(BorderFactory.createEmptyBorder(30, //top
		30, //left
		60, //bottom
		30) //right
		);
		pane.setLayout(new GridLayout(0, 1));
		pane.add(button1);
		pane.add(button2);
		pane.add(button3);
		pane.add(button4);
		pane.add(button5);
		pane.add(button6);
		pane.add(button7);
		pane.add(button8);
		pane.add(button9);

		return pane;
	}

	public void createFrame(ApplicationComponent app) {
		_app = app;
		//Create the top-level container and add contents to it.
		_frm = new JFrame("SwingApplication");
		Visualizer2DControler vc = new Visualizer2DControler();
		Component contents = vc.createComponents();
		_frm.getContentPane().add(contents, BorderLayout.CENTER);

		//Finish setting up the frame, and show it.
		_frm.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				//do nothing, just close the window 
			}
		});
		_frm.pack();
		_frm.setVisible(true);
	}
}

`timescale 1ns / 1ps

////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer:
//
// Create Date:   21:38:07 05/07/2015
// Design Name:   test
// Module Name:   C:/Users/Harsh/Desktop/VSD/prev_scr/test_tb.v
// Project Name:  prev_scr
// Target Device:  
// Tool versions:  
// Description: 
//
// Verilog Test Fixture created by ISE for module: test
//
// Dependencies:
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
////////////////////////////////////////////////////////////////////////////////

module test_tb;

	// Inputs
	reg clk;
	reg [13:0] R_fixed;
	reg [13:0] angle1;
	reg [3:0] adder;
	reg [4:0] error;

	// Outputs
	wire [3:0] X;
	wire [3:0] Y;

	// Instantiate the Unit Under Test (UUT)
	test uut (
		.clk(clk), 
		.R_fixed(R_fixed), 
		.angle1(angle1), 
		.adder(adder), 
		.X(X), 
		.Y(Y), 
		.error(error)
	);

	always @(clk)
	begin
		#0.001 clk <= ~clk;
	end

	initial begin
		// Initialize Inputs
		clk = 0;
		R_fixed = 14'b00011110000000 ;
		angle1 = 14'b00011100000000;
		adder = 4'b0100;
		error = 5'b00111;

		// Wait 100 ns for global reset to finish
		#100;
        
		// Add stimulus here

	end
      
endmodule


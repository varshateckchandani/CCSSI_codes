`timescale 1ns / 1ps

////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer:
//
// Create Date:   20:57:22 05/07/2015
// Design Name:   main
// Module Name:   C:/Users/Harsh/Desktop/VSD/CCSSI_final_group2/main_tb.v
// Project Name:  CCSSI_final_group2
// Target Device:  
// Tool versions:  
// Description: 
//
// Verilog Test Fixture created by ISE for module: main
//
// Dependencies:
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
////////////////////////////////////////////////////////////////////////////////

module main_tb;

	// Inputs
	reg clk;
	reg [13:0] R_fixed;
	reg [13:0] angle1;
	reg [3:0] adder;
	reg [4:0] error;
	reg [13:0] angle2;

	// Outputs
	wire [7:0] final_cord1;
	wire [7:0] final_cord2;

	// Instantiate the Unit Under Test (UUT)
	main uut (
		.clk(clk), 
		.R_fixed(R_fixed), 
		.angle1(angle1), 
		.adder(adder), 
		.error(error), 
		.angle2(angle2), 
		.final_cord1(final_cord1), 
		.final_cord2(final_cord2)
	);
	always @(clk)
	begin
		#0.001 clk <= ~clk;
	end


	initial begin
		// Initialize Inputs
		clk = 0;
		R_fixed = 14'b00011111100000;
		angle1 = 14'b01001100000000;
		angle2 = 14'b01101000000000;
		adder = 4'b1000;
		error=5'b00111;
		
		// Wait 100 ns for global reset to finish
		#100;
        
		// Add stimulus here

	end
      
endmodule


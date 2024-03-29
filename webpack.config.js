const path = require('path');

module.exports = {
	entry: './src/index.ts',
	devtool: 'inline-source-map',
	mode: 'development',
	module: {
		rules: [
			{
				test: /\.tsx?$/,
				use: [
					'ts-loader'
				],
				exclude: /node_modules/,
			},
		],
	},
	resolve: {
		extensions: ['.tsx', '.ts', '.js'],

		fallback: {
			"crypto": false,
			"path": false,
			"fs": false
		}
	},
	output: {
		filename: 'bundle.js',
		path: path.resolve(__dirname, 'basic-arap'),
	},
	watch: true
};
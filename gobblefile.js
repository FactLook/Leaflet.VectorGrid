

var gobble = require('gobble');

var concatenatedJs = gobble('src').transform('concat', {
	dest: 'Leaflet.VectorGrid.js',
	files: [
		'Leaflet.Renderer.SVG.Tile.js',
		'Leaflet.Renderer.Canvas.Tile.js',
		'Leaflet.VectorGrid.js',
		'Leaflet.VectorGrid.Slicer.js',
		'Leaflet.VectorGrid.Protobuf.js'
	]
});

var vendor = gobble('vendor');

var bundled = gobble([ vendor, concatenatedJs ]).transform('concat', {
	dest: 'Leaflet.VectorGrid.bundled.js',
	files: [
		'geojson-vt-dev.js',
		'pbf-dev.js',
		'topojson.js',
		'vector-tile-dev.js',
		'Leaflet.VectorGrid.js'
	],
// 	writeSourcemap: false
});

var demo = gobble('demo').moveTo('demo');
var leaflet = vendor.include([
	'leaflet-src.js',
	'leaflet-src.map',
	'leaflet.css'
]).moveTo('demo');


module.exports = gobble([
	concatenatedJs,
	concatenatedJs.transform('uglifyjs', { ext: '.min.js' }),
	bundled,
	bundled.transform('uglifyjs', { ext: '.min.js' }),
	demo,
	leaflet
]);

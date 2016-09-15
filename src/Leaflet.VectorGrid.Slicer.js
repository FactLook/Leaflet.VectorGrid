


// geojson-vt powered!
// NOTE: Assumes the global `geojsonvt` exists!!!
L.VectorGrid.Slicer = L.VectorGrid.extend({

	options: {
		vectorTileLayerName: 'sliced',
		extent: 4096,	// Default for geojson-vt
		maxZoom: 14  	// Default for geojson-vt
	},

	initialize: function(geojson, options) {
		L.VectorGrid.prototype.initialize.call(this, options);


			this._slicers = {};
        this.addJsonLayer(geojson);
		/*if (geojson.type && geojson.type === 'Topology') {
			// geojson is really a topojson
			for (var layerName in geojson.objects) {
				this._slicers[layerName] = geojsonvt(
					topojson.feature(geojson, geojson.objects[layerName])
				, this.options);
                console.log('topojson layer:', layerName);
			}
		} else {
			// For a geojson, create just one vectortilelayer named with the value
			// of the option.
			// Inherits available options from geojson-vt!
                console.log('Geojson layer:');
			this._slicers[this.options.vectorTileLayerName] = geojsonvt(geojson, this.options);
		}*/
	},

  addJsonLayer : function(geojson){
        if (geojson.type && geojson.type === 'Topology') {
            // geojson is really a topojson
            for (var layerName in geojson.objects) {
                var existingLayer = this._slicers[layerName];
                if(!existingLayer){
                    this._slicers[layerName] = geojsonvt(
                        topojson.feature(geojson, geojson.objects[layerName])
                        , this.options);
                    //console.log('topojson layer:', layerName);
                }
            }
        } else {
            // For a geojson, create just one vectortilelayer named with the value
            // of the option.
            // Inherits available options from geojson-vt!
            //console.log('Geojson layer:');
            var existingLayer = this._slicers[layerName];
            if(!existingLayer){
                this._slicers[this.options.vectorTileLayerName] = geojsonvt(geojson, this.options);
            }
        }
    },

	_getVectorTilePromise: function(coords) {

		var tileLayers = {};

		for (var layerName in this._slicers) {
			var slicer = this._slicers[layerName];
			var slicedTileLayer = slicer.getTile(coords.z, coords.x, coords.y);

// 			console.log(coords, slicedTileLayer && slicedTileLayer.features && slicedTileLayer.features.length || 0);

			if (slicedTileLayer) {
				var vectorTileLayer = {
					features: [],
					extent: this.options.extent,
					name: this.options.vectorTileLayerName,
					length: slicedTileLayer.features.length
				}

				for (var i in slicedTileLayer.features) {
					var feat = {
						geometry: slicedTileLayer.features[i].geometry,
						properties: slicedTileLayer.features[i].tags,
						type: slicedTileLayer.features[i].type	// 1 = point, 2 = line, 3 = polygon
					}
					vectorTileLayer.features.push(feat);
				}

				tileLayers[layerName] = vectorTileLayer;
			}

		}

		return new Promise(function(resolve){ return resolve({ layers: tileLayers })});
	},

});


L.vectorGrid.slicer = function (geojson, options) {
	return new L.VectorGrid.Slicer(geojson, options);
};


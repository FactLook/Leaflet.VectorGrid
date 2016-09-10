L.VectorGrid = L.GridLayer.extend({

	options: {
		rendererFactory: L.svg.tile,
		vectorTileLayerStyles: {}
	},

	createTile: function(coords, done) {
		var renderer = this.options.rendererFactory(this.getTileSize(), this.options);

		var vectorTilePromise = this._getVectorTilePromise(coords);

		vectorTilePromise.then( function renderTile(vectorTile) {

			for (var layerName in vectorTile.layers) {
				var layer = vectorTile.layers[layerName];

				/// NOTE: THIS ASSUMES SQUARE TILES!!!!!1!
				var pxPerExtent = this.getTileSize().x / layer.extent;

        var layerStyle;
        if (this.options.vectorTileLayerStyles instanceof Function) {
          layerSTyle = this.options.vectorTileLayerStyles;
        } else if (this.options.vectorTileLayerStyles instanceof Object) {
          layerStyle = this.options.vectorTileLayerStyles[layerName]
        } else {
          layerStyle = L.Path.prototype.options;
        }

        var onClick = getEventHandler(layerName, this.options.onClick);
        var onMouseDown = getEventHandler(layerName, this.options.onMouseDown);
        var onMouseUp = getEventHandler(layerName, this.options.onMouseUp);
        var onMouseOver = getEventHandler(layerName, this.options.onMouseOver);
        var onMouseMove = getEventHandler(layerName, this.options.onMouseMove);
        var onMouseOut = getEventHandler(layerName, this.options.onMouseOut);

				for (var i in layer.features) {
					var feat = layer.features[i];

					if (feat.type > 1) { // Lines, polygons

						this._mkFeatureParts(feat, pxPerExtent);

						/// Style can be a callback that is passed the feature's
						/// properties and tile zoom level...
						var styleOptions = (layerStyle instanceof Function) ?
						layerStyle(feat.properties, coords.z) :
						layerStyle;

						if (!(styleOptions instanceof Array)) {
							styleOptions = [styleOptions];
						}

            feat._map = this._map;
						feat.getLayerName = function() { return layerName; };
						feat.getCoords = function() { return coords; };

            feat._onClick = onClick;
            feat._onMouseMove = onMouseMove;
            feat._onMouseDown = onMouseDown;
            feat._onMouseUp = onMouseUp;
            feat._onMouseOver = onMouseOver;
            feat._onMouseOut = onMouseOut;

						/// Style can be an array of styles, for styling a feature
						/// more than once...
						for (var j in styleOptions) {
							var style = L.extend({}, L.Path.prototype.options, styleOptions[j]);

							if (feat.type === 1) { // Points
								style.fill = false;
							} else if (feat.type === 2) {	// Polyline
								style.fill = false;
							}

							feat.options = style;
							renderer._initPath( feat );
							renderer._updateStyle( feat );

							if (feat.type === 1) { // Points
								// 							style.fill = false;
							} else if (feat.type === 2) {	// Polyline
								style.fill = false;
								renderer._updatePoly(feat, false);
							} else if (feat.type === 3) {	// Polygon
								renderer._updatePoly(feat, true);
							}

							renderer._addPath( feat );
						}

					} else {
						// Feat is a point (type === 1)

						/// FIXME!!!
					}
				}

			}
			L.Util.requestAnimFrame(done);
		}.bind(this));

		return renderer.getContainer();
	},

	// Fills up feat._parts based on the geometry and pxPerExtent,
	// pretty much as L.Polyline._projectLatLngs and L.Polyline._clipPoints
	// would do but simplified as the vectors are already simplified/clipped.
	_mkFeatureParts: function(feat, pxPerExtent) {

		var rings = feat.geometry;

		feat._parts = [];
		for (var i in rings) {
			var ring = rings[i];
			var part = [];
			for (var j in ring) {
				var coord = ring[j];
				if ('x' in coord) {
					// Protobuf vector tiles return {x: , y:}
					part.push(L.point(coord.x * pxPerExtent, coord.y * pxPerExtent));
				} else {
					// Geojson-vt returns [,]
					part.push(L.point(coord[0] * pxPerExtent, coord[1] * pxPerExtent));
				}
			}
			feat._parts.push(part);
		}

	},

});

L.vectorGrid = function (options) {
	return new L.VectorGrid(options);
};

function getEventHandler(layerName, handler) {
  if (typeof handler === 'function') {
    return handler;
  } else if (typeof handler === 'object') {
    if (typeof handler[layerName] === 'function') {
      return handler[layerName];
    }
  }
}



L.SVG.Tile = L.SVG.extend({

	initialize: function (tileSize, options) {
		L.SVG.prototype.initialize.call(this, options);
		this._size = tileSize;

		this._initContainer();
		this._container.setAttribute('width', this._size.x);
		this._container.setAttribute('height', this._size.y);
		this._container.setAttribute('viewBox', [0, 0, this._size.x, this._size.y].join(' '));
	},

	getContainer: function() {
		return this._container;
	},

// 	onAdd: function() {},
	onAdd: L.Util.FalseFn,

	_initContainer: function() {
		L.SVG.prototype._initContainer.call(this);
		var rect =  L.SVG.create('rect');

// 		rect.setAttribute('x', 0);
// 		rect.setAttribute('y', 0);
// 		rect.setAttribute('width', this._size.x);
// 		rect.setAttribute('height', this._size.y);
// 		rect.setAttribute('fill', 'transparent');
// 		rect.setAttribute('stroke', 'black');
// 		rect.setAttribute('stroke-width', 2);
// 		this._rootGroup.appendChild(rect);
	},

	/// TODO: Modify _initPath to include an extra parameter, a group name
	/// to order symbolizers by z-index
	_initPath: function(layer) {
		L.SVG.prototype._initPath.call(this, layer);

		var path = layer._path;

		if (typeof layer.onClick === 'function') {
			if (layer.type === 2) {
				path.setAttribute('pointer-events', 'stroke');
			} else if (layer.type === 3) {
				path.setAttribute('pointer-events', 'fill');
			}

			path.addEventListener('click', function(e) {
				layer.onClick({ type: 'click', target: layer });
			});
		}
	},

	_addPath: function (layer) {
		this._rootGroup.appendChild(layer._path);
	},

});


L.svg.tile = function(tileSize, opts){
	return new L.SVG.Tile(tileSize, opts);
}

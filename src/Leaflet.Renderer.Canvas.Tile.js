

L.Canvas.Tile = L.Canvas.extend({

	initialize: function (tileSize, options) {
		L.Canvas.prototype.initialize.call(this, options);
		this._size = tileSize;

		this._initContainer();
		this._container.setAttribute('width', this._size.x);
		this._container.setAttribute('height', this._size.y);
		this._layers = {};
		this._drawnLayers = {};
	},

	getContainer: function() {
		return this._container;
	},

	onAdd: L.Util.FalseFn,

	_initContainer: function () {
		var container = this._container = document.createElement('canvas');

// 		L.DomEvent
// 			.on(container, 'mousemove', L.Util.throttle(this._onMouseMove, 32, this), this)
// 			.on(container, 'click dblclick mousedown mouseup contextmenu', this._onClick, this)
// 			.on(container, 'mouseout', this._handleMouseOut, this);

		this._ctx = container.getContext('2d');
	},

	/// TODO: Modify _initPath to include an extra parameter, a group name
	/// to order symbolizers by z-index
	_initPath: function(layer) {
		L.Canvas.prototype._initPath.call(this, layer);

		var path = layer._path;

		if (typeof layer.onClick === 'function') {
			if (layer.type === 2) {
				path.setAttribute('pointer-events', 'stroke');
			} else if (layer.type === 3) {
				path.setAttribute('pointer-events', 'fill');
			} else if (layer.type === 1) {
        path.setAttribute('pointer-events', 'painted');
      }

			path.addEventListener('click', function(e) {
				layer.onClick({ type: 'click', target: layer });
			});
		}
	},
});

L.canvas.tile = function(tileSize, opts){
	return new L.Canvas.Tile(tileSize, opts);
}

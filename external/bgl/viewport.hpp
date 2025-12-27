#pragma once

namespace bgl {

	class Viewport {
	public:
		inline constexpr Viewport(int width, int height) : m_width(width), m_height(height) {}

		inline constexpr int width() const { return m_width; }

		inline constexpr int height() const { return m_height; }

		template<bool CLAMPED>
		inline constexpr void clip_space_to_screen_inv_w(float& x, float& y, float& w) const {

			w = 1.0f / w;

			x *= w;
			y *= w;

			x = x * 0.5f + 0.5f;
			y = y * 0.5f + 0.5f;

			if constexpr (CLAMPED) {
				if (x <= 0.f)
					x = 0.f;
				if (x >= 1.f)
					x = 1.f;

				if (y <= 0.f)
					y = 0.f;
				if (y >= 1.f)
					y = 1.f;
			}

			x = x * width();
			y = y * height();
		}

		template<bool CLAMPED>
		inline constexpr void clip_space_to_screen(float& x, float& y, float w) const {
			clip_space_to_screen_inv_w<CLAMPED>(x, y, w);
		}

		inline constexpr void clip_space_to_screen(float& x, float& y, float w) const {
			clip_space_to_screen_inv_w<true>(x, y, w);
		}

		inline constexpr void clip_space_to_screen_unclamped(float& x, float& y, float w) const {
			clip_space_to_screen_inv_w<false>(x, y, w);
		}

		inline constexpr void clip_space_to_screen_inv_w(float& x, float& y, float& w) const {
			clip_space_to_screen_inv_w<true>(x, y, w);
		}


		inline constexpr void clip_space_to_screen_inv_w_unclamped(float& x, float& y, float& w) const {
			clip_space_to_screen_inv_w<false>(x, y, w);
		}

	private:

		int m_width;
		int m_height;
	};

}